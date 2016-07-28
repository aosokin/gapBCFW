function [model, gap_vec, frac_pass, cache, data_pass_info] = onePass_BCFW_hybrid(param, options, model, gap_vec, cache, gap_value_for_cache)
% [model, gap_vec, frac_pass, cache, data_pass_info] = onePass_BCFW_hybrid(param, options, model, gap_vec, cache, gap_value_for_cache)
%
% onePass_BCFW_hybrid does one pass over the dataset (samples n objects) for all the variants of BCFW algorithm.
% Designed to be called from solver_BCFW_hybrid.m

%% get problem description
phi = param.featureFn; % for \phi(x,y) feature mapping
loss = param.lossFn; % for L(ytruth, ypred) loss function
oracle = param.oracleFn;
n = param.n; % number of training examples
lambda = options.lambda;

if ~exist('model', 'var')
    model = initialize_model_zeros(param);
end
if ~exist('gap_vec', 'var')
    gap_vec = inf(param.n,1);
end
if ~exist('gap_value_for_cache', 'var')
    gap_value_for_cache = inf;
end

% do removeFromActiveSetIfFullCache only if there no away steps and no pairwise steps
removeFromActiveSetIfFullCache = ~options.use_away_steps && ~options.use_pairwise_steps;

%% check the positivity constraints
if isfield(param, 'positivity')
    positivity = param.positivity;
    % check that the constraint are respected
    if ~all(model.w(param.positivity==1)>=0)
        fprintf('CAUTION (%s): The model.w does not respect %d positivity constraints\n', mfilename, sum(model.w(param.positivity==1)<0));
    end
    % need to reconstruct the non truncated version of v from the model
    % keep a non trucanted version of w (A\alpha = v)
    model.v = sum(model.wMat, 2);
    % check if sparsity (note that v (as w) is full).
    if issparse(model.v)
        model.v = full(model.v);
    end
    % making sure that the quantities are consistent
    v_truncated = model.v;
    v_truncated(param.positivity==1) = subplus(model.v(param.positivity==1));
    % check if the reconstruction of v made sense
    if norm(model.w - v_truncated) >= 1e-12
        fprintf('CAUTION (%s): The model.w does not match the truncated sum of model.wMat. Error: %f\n', mfilename, norm(model.w - v_truncated));
        model.w = v_truncated;
    end 
    clear v_truncated;
else
    positivity = [];
end

%% Initialize the pass
frac_pass = 0;
perm = [];
if (isequal(options.sample, 'perm'))
    perm = randperm(n);
end

data_pass_info = struct;
data_pass_info.numCacheHitsThisDataPass = 0;
data_pass_info.noCacheHit = 0;
data_pass_info.globalGapHit = 0;
data_pass_info.localGapHit = 0;
data_pass_info.doubleGapHit = 0;
data_pass_info.numFwSteps = 0;
data_pass_info.numAwaySteps = 0;
data_pass_info.numPairwiseSteps = 0;

number_drop_steps = zeros(n,1); % record the number of drop steps for each block

%% One effective pass: visiting n objects 
for dummy = 1 : n
    
    % exit the loop if heuristic gap is small enough
    if options.quit_passes_heuristic_gap && (sum(gap_vec) <= options.gap_threshold*options.quit_passes_heuristic_gap_eps_multiplyer)
        break;
    end
    
    %% Picking random example:
    if (isequal(options.sample, 'uniform'))
        i = randi(n); % uniform sampling
    elseif(isequal(options.sample, 'gap'))
        if(min(gap_vec)<-1E-5)
            fprintf('CAUTION (%s): gap is negative (min(gaps)):%.6f\n', mfilename, min(gap_vec));
        end
        gap_vec(gap_vec<0) = 0;
        if ~isfield(options, 'gapDegree')
            probs = gap_vec;
        else
            probs = gap_vec .^ options.gapDegree;
        end
        i = randsample_fast(probs); % gap based sampling
    elseif(isequal(options.sample, 'maxGap'))
        if( all(isinf(gap_vec)) )
            i = randi(n);
        else
            if(min(gap_vec)<-1E-5)
                fprintf('CAUTION (%s): gap is negative (min(gaps)):%.6f\n', mfilename, min(gap_vec));
            end
            gap_vec(gap_vec<0) = 0;
            [~, i] = max(gap_vec); % gap based sampling
        end
    elseif (isequal(options.sample, 'perm'))
        i = perm(dummy); % random permutation
    else
        error([mfilename, 'Unknown gap sampling method!']);
    end
    
    %% Get the FW corner
    if options.useCache
        % try to solve the loss-augmented inference for point i using cache
        
        loss_augmented_scores = cache_get_scores( cache{i}, model.w );
        [~, cache_best_index] = max(loss_augmented_scores);
        
        [ psi_cache, loss_cache, yhash_i ] = cache_get_entry_by_index( cache{i}, cache_best_index );
        
        w_cache = 1/(n*lambda) * psi_cache;
        ell_cache = (1/n)* loss_cache;
        cache_maximal_improvement = lambda*(model.w'*(model.wMat(:,i) - w_cache)) - (model.ellMat(i) - ell_cache);
        
        if ~removeFromActiveSetIfFullCache
            % if throwing away point from active set is not allowed this should hold
            %assert(gap_cache > -100*eps);
            if (cache_maximal_improvement < -100*eps)
                fprintf('CAUTION: gap-cache assertion of %s fails. Gap value = %f\n', mfilename, cache_maximal_improvement);
            end
        end
        
        % criterion for hitting the cache
        cache_hit_threshold = max( gap_vec(i) * options.cacheFactor, options.cacheNu * gap_value_for_cache/n);
        cache_hit_flag = cache_maximal_improvement >= cache_hit_threshold;
        
        if ~cache_hit_flag
            data_pass_info.noCacheHit = data_pass_info.noCacheHit + 1;
        else
            if gap_vec(i) * options.cacheFactor >= cache_hit_threshold && options.cacheNu * gap_value_for_cache/n >= cache_hit_threshold
                data_pass_info.doubleGapHit = data_pass_info.doubleGapHit + 1;
            elseif gap_vec(i) * options.cacheFactor >= cache_hit_threshold
                data_pass_info.localGapHit = data_pass_info.localGapHit + 1;
            elseif options.cacheNu * gap_value_for_cache/n >= cache_hit_threshold
                data_pass_info.globalGapHit = data_pass_info.globalGapHit + 1;
            else
                error('Something weird has happened with the cache hit!');
            end
        end
        
    else
        cache_maximal_improvement = -inf;
        cache_hit_flag = false;
    end
    if ~options.useCache || ~cache_hit_flag % if cache point is not good enough
        % solve the loss-augmented inference for point i
        ystar_i = oracle(param, model, param.patterns{i}, param.labels{i});
        frac_pass = frac_pass + 1/n;
        
        % define the update quantities:
        % [note that lambda*w_s is subgradient of 1/n*H_i(w) ]
        % psi_i(y) := phi(x_i,y_i) - phi(x_i, y)
        psi_i =   phi(param, param.patterns{i}, param.labels{i}) ...
            - phi(param, param.patterns{i}, ystar_i);
        if isfield(param, 'cache_type_single') && param.cache_type_single
            % if cache is going to be single we need to truncate to single precision to avoid future in
            psi_i = single(psi_i);
            % need to cast back to double to have all the computations in double
            psi_i = double(psi_i);
        end
        w_s = 1/(n*lambda) * psi_i;
        loss_i = loss(param, param.labels{i}, ystar_i);
        ell_s = 1/n*loss_i;
        gap_i_FW = lambda*(model.w'*(model.wMat(:,i) - w_s)) - (model.ellMat(i) - ell_s);
        
        % sanity check, if this assertion fails, probably there is a bug in the
        % maxOracle or in the featuremap
        if ((loss_i - model.w'*psi_i) < -1e-12)
            fprintf('CAUTION: oracles assertion of %s fails on %f\n', mfilename, loss_i - model.w'*psi_i);
            ystar_i = param.labels{i};
            psi_i = zeros(size(psi_i));
            w_s = 1/(n*lambda) * psi_i;
            loss_i = loss(param, param.labels{i}, ystar_i);
            ell_s = 1/n*loss_i;
            gap_i_FW = lambda*(model.w'*(model.wMat(:,i) - w_s)) - (model.ellMat(i) - ell_s);
        end
        
        % we update the gap for this example; used for gap sampling, need this to be non-negative
        gap_vec(i) = max(gap_i_FW,0);
        
        maximal_improvement_FW = gap_i_FW;
        if options.update_dual_vars
            yhash_i = param.hashFn(ystar_i);
        end
    else % cache point is good enough for update
        % define the update quantities using the cache point
        % [note that lambda*w_s is subgradient of 1/n*H_i(w) ]
        % psi_i(y) := phi(x_i,y_i) - phi(x_i, y)
        
        psi_i = psi_cache;
        w_s = w_cache;
        loss_i = loss_cache;
        ell_s = ell_cache;
        maximal_improvement_FW = cache_maximal_improvement;
        
        data_pass_info.numCacheHitsThisDataPass = data_pass_info.numCacheHitsThisDataPass + 1;
    end
    
    %% Get the away corner
    if options.use_away_steps || options.use_pairwise_steps
        alphas = cache_get_dual_vars( cache{i} );
        active_set_indices = find( alphas > 0 );
        if ~options.useCache
            active_set_scores = cache_get_scores( cache{i}, model.w, active_set_indices );
            % update the loss_augmented_scores if they were not computed before
            loss_augmented_scores = nan(numel(alphas), 1);
            loss_augmented_scores(active_set_indices) = active_set_scores;
        else
            % loss_augmented_score is already computed for the cache
            active_set_scores = loss_augmented_scores(active_set_indices);
        end
        [~, cache_away_index] = min(active_set_scores);
        cache_away_index = active_set_indices(cache_away_index);
        
        [ psi_away, loss_away, yhash_away ] = cache_get_entry_by_index( cache{i}, cache_away_index );
        w_a = 1/(n*lambda) * psi_away;
        ell_a = (1/n)* loss_away;
        alpha_a = alphas(cache_away_index);
        
        gap_away = lambda*(model.w'*(w_a-model.wMat(:,i))) + model.ellMat(i) - ell_a;
        if gap_away <= -100*eps
            fprintf('CAUTION: gap_away assetion in %s fails on %f\n', mfilename, gap_away);
        end    
    end
    
    %% the current objective function
    dual_value_current = -objective_function( model.w, model.ell, lambda);
    
    %% try the FW step
    if options.use_FW_steps
        % line_search:
        gamma_opt_FW_step = maximal_improvement_FW / (lambda*( (model.wMat(:,i) - w_s)'*(model.wMat(:,i) - w_s) +eps));
        
        % truncation to [0,1] with the update of dual variables
        if gamma_opt_FW_step >= 1
            is_drop_step_FW = 1;
            gamma_FW_step = 1;
            number_drop_steps(i) = number_drop_steps(i) + 1;
        else
            gamma_FW_step = max(0, gamma_opt_FW_step);
            is_drop_step_FW = 0;
        end
        % try to do the FW step
        if gamma_FW_step > 0
            % update the weights and ell variables
            if ~isempty(positivity)
                v_FW = model.v - model.wMat(:,i); % this is v^(k)-w_i^(k)
                wMat_i_FW = (1-gamma_FW_step)*model.wMat(:,i) + gamma_FW_step*w_s;
                v_FW = v_FW + wMat_i_FW; % this is v^(k+1) = v^(k)-w_i^(k)+w_i^(k+1)
                % check for sparsity (sanity check)
                if issparse(v_FW)
                    v_FW = full(v_FW);
                end
                % truncate w (optimize over \beta)
                w_FW = v_FW;
                w_FW(positivity==1) = subplus(v_FW(positivity==1));
            else
                w_FW = model.w - model.wMat(:,i); % this is w^(k)-w_i^(k)
                wMat_i_FW = (1-gamma_FW_step)*model.wMat(:,i) + gamma_FW_step*w_s;
                w_FW = w_FW + wMat_i_FW; % this is w^(k+1) = w^(k)-w_i^(k)+w_i^(k+1)
                %check for sparsity
                if issparse(w_FW)
                    w_FW = full(w_FW);
                end
            end
            ell_FW = model.ell - model.ellMat(i); % this is ell^(k)-ell_i^(k)
            ellMat_i_FW = (1-gamma_FW_step)*model.ellMat(i) + gamma_FW_step*ell_s;
            ell_FW = ell_FW + ellMat_i_FW; % this is ell^(k+1) = ell^(k)-ell_i^(k)+ell_i^(k+1)
            
            dual_value_FW = -objective_function( w_FW, ell_FW, lambda);
            improvement_FW_step = dual_value_FW - dual_value_current;
            
            if improvement_FW_step < -eps
                fprintf('CAUTION: FW step in %s leads to negative improvement: %f! Something is wrong with the FW step.', mfilename, improvement_FW_step);
            end
        else
            improvement_FW_step = 0;
        end
    else
        improvement_FW_step = -inf;
    end
    
    %% try the pairwise step
    if options.use_pairwise_steps
        max_step_size_pairwise_step = alpha_a;
        % line_search:
        gamma_opt_pairwise_step = (gap_away+maximal_improvement_FW) / (lambda * norm(w_s - w_a)^2 + eps);
        if gamma_opt_pairwise_step >= max_step_size_pairwise_step
            % we have a drop step!
            gamma_pairwise_step = max_step_size_pairwise_step;
            number_drop_steps(i) = number_drop_steps(i) + 1;
            is_drop_step_pairwise = 1;
        else
            gamma_pairwise_step = max(0, gamma_opt_pairwise_step);
            is_drop_step_pairwise = 0;
        end
        % try to do the FW step
        if gamma_pairwise_step > 0
            % update the weights and ell variables
            if ~isempty(positivity)
                v_pairwise = model.v - model.wMat(:,i); % this is v^(k)-w_i^(k)
                wMat_i_pairwise = model.wMat(:,i) + gamma_pairwise_step*(w_s-w_a);
                v_pairwise = v_pairwise + wMat_i_pairwise; % this is v^(k+1) = v^(k)-w_i^(k)+w_i^(k+1)
                
                if issparse(v_pairwise)
                    v_pairwise = full(v_pairwise);
                end
                
                % truncate w (optimize over \beta)
                w_pairwise = v_pairwise;
                w_pairwise(positivity==1) = subplus(v_pairwise(positivity==1));
            else
                w_pairwise = model.w - model.wMat(:,i); % this is w^(k)-w_i^(k)
                wMat_i_pairwise = model.wMat(:,i) + gamma_pairwise_step*(w_s-w_a);
                w_pairwise = w_pairwise + wMat_i_pairwise; % this is w^(k+1) = w^(k)-w_i^(k)+w_i^(k+1)
                if issparse(w_pairwise)
                    w_pairwise = full(w_pairwise);
                end
            end
            ell_pairwise = model.ell - model.ellMat(i); % this is ell^(k)-ell_i^(k)
            ellMat_i_pairwise = model.ellMat(i) + gamma_pairwise_step*(ell_s-ell_a);
            ell_pairwise = ell_pairwise + ellMat_i_pairwise; % this is ell^(k+1) = ell^(k)-ell_i^(k)+ell_i^(k+1)
            
            dual_value_pairwise = -objective_function( w_pairwise, ell_pairwise, lambda);
            improvement_pairwise_step = dual_value_pairwise - dual_value_current;
            
            if improvement_pairwise_step < -eps
                fprintf('CAUTION: Pairwise step in %s leads to negative improvement: %f! Something is wrong with the pairwise step.', mfilename, improvement_pairwise_step);
            end
        else
            improvement_pairwise_step = 0;
        end
    else
        improvement_pairwise_step = -inf;
    end
    
    %% try the away step
    if options.use_away_steps
        % only do away step if active set has more than 1 active elements:
        number_active = sum( alphas > 0 );
        if number_active > 1
            % can consider doing an away step
            max_step_size_away_step = alpha_a / (1-alpha_a);
            is_drop_step_away = 0;
            
            % line_search:
            gamma_opt_away_step = gap_away / (lambda * norm(model.wMat(:,i) - w_a)^2 + eps);
            if gamma_opt_away_step >= max_step_size_away_step
                % we have a drop step!
                gamma_away_step = max_step_size_away_step;
                number_drop_steps(i) = number_drop_steps(i) + 1;
                is_drop_step_away = 1;
            else
                gamma_away_step = max(0, gamma_opt_away_step);
            end
            
            if gamma_away_step > 0
                % update the weights and ell variables
                if ~isempty(positivity)
                    v_away = model.v - model.wMat(:,i); % this is v^(k)-w_i^(k)
                    wMat_i_away = (1+gamma_away_step)*model.wMat(:,i) - gamma_away_step*w_a;
                    v_away = v_away + wMat_i_away; % this is v^(k+1) = v^(k)-w_i^(k)+w_i^(k+1)
                    % check for sparsity
                    if issparse(v_away)
                        v_away = full(v_away);
                    end
                    % truncate w (optimize over \beta)
                    w_away = v_away;
                    w_away(positivity==1) = subplus(v_away(positivity==1));
                else
                    w_away = model.w - model.wMat(:,i); % this is w^(k)-w_i^(k)
                    wMat_i_away = (1+gamma_away_step)*model.wMat(:,i) - gamma_away_step*w_a;
                    w_away = w_away + wMat_i_away; % this is w^(k+1) = w^(k)-w_i^(k)+w_i^(k+1)
                    % check for sparsity
                    if issparse(w_away)
                        w_away = full(w_away);
                    end
                end
                ell_away = model.ell - model.ellMat(i); % this is ell^(k)-ell_i^(k)
                ellMat_i_away = (1+gamma_away_step)*model.ellMat(i) - gamma_away_step*ell_a;
                ell_away = ell_away + ellMat_i_away; % this is ell^(k+1) = ell^(k)-ell_i^(k)+ell_i^(k+1)
                
                dual_value_away = -objective_function( w_away, ell_away, lambda);
                improvement_away_step = dual_value_away - dual_value_current;
                
                if improvement_away_step < -eps
                    fprintf('CAUTION: Away step in %s leads to negative improvement: %f! Something is wrong with the away step.', mfilename, improvement_away_step);
                end
            else
                improvement_away_step = 0;
            end
        else
            % can't try away step
            improvement_away_step = -inf;
        end
    else
        improvement_away_step = -inf;
    end
    
    %% choose the best step in a greedy way
    best_step_type = '';
    best_step_improvement = max( [improvement_FW_step; improvement_pairwise_step; improvement_away_step] );
    if improvement_away_step >= best_step_improvement
        best_step_type = 'away';
        is_drop_step = is_drop_step_away;
    end
    if improvement_pairwise_step >= best_step_improvement
        best_step_type = 'pairwise';
        is_drop_step = is_drop_step_pairwise;
    end
    if improvement_FW_step >= best_step_improvement
        best_step_type = 'FW';
        is_drop_step = is_drop_step_FW;
    end
    
    %% do the primal step
    % only if it leads to a non-zero improvement
    if best_step_improvement > 0
        switch best_step_type
            case 'FW'
            model.w = w_FW;
            model.ell = ell_FW;
            model.wMat(:,i) = wMat_i_FW;
            model.ellMat(i) = ellMat_i_FW;
            if ~isempty(positivity)
                model.v = v_FW;
            end
            data_pass_info.numFwSteps = data_pass_info.numFwSteps + 1;
        case 'away'
            model.w = w_away;
            model.ell = ell_away;
            model.wMat(:,i) = wMat_i_away;
            model.ellMat(i) = ellMat_i_away;
            if ~isempty(positivity)
                model.v = v_away;
            end
            data_pass_info.numAwaySteps = data_pass_info.numAwaySteps + 1;
        case 'pairwise'
            model.w = w_pairwise;
            model.ell = ell_pairwise;
            model.wMat(:,i) = wMat_i_pairwise;
            model.ellMat(i) = ellMat_i_pairwise;
            if ~isempty(positivity)
                model.v = v_pairwise;
            end
            data_pass_info.numPairwiseSteps = data_pass_info.numPairwiseSteps + 1;
        otherwise
            error('ERROR: internal error: unknown type of step!')
        end
    end
    
    %% even if there is no improvement save the variable, in case of an oracle call
    if options.update_dual_vars
        % try to get object from cache, if it is not there it will be added
        % need the two calls here to make sure the two objects are in the cache
        if options.use_away_steps || options.use_pairwise_steps
            [cache{i}, cache_away_index ] = cache_get_entry_by_hash( cache{i}, yhash_away, psi_away, loss_away, ...
                loss_augmented_scores, options.maxCacheSize, removeFromActiveSetIfFullCache );
        end
        [cache{i}, cache_towards_index ] = cache_get_entry_by_hash( cache{i}, yhash_i, psi_i, loss_i, ...
            loss_augmented_scores, options.maxCacheSize, removeFromActiveSetIfFullCache );
        
        if (options.use_away_steps || options.use_pairwise_steps) 
            % check that FW and away corners are different
            % - this condition might break if active set is of size one
            % - this condition might also break if active_set_scores are exactly equal for multiple elements
            %   this can happen when the point was already optimal
            if (cache_towards_index == cache_away_index) && (numel(active_set_indices) > 1) && (improvement_pairwise_step > 0 || improvement_away_step > 0)
                fprintf('CAUTION! Away and FW corners are the same in %s. It is either a bug or some strange edge case: might be worth investigating.\n', mfilename);
            end
        end
    end
    
    %% do the best step
    if options.update_dual_vars
        if best_step_improvement > 0
            % update the alpha variables 
            alphas = cache_get_dual_vars( cache{i} ); % need to get the alphas again, because cache updates might have changed them
            switch best_step_type
                case 'FW'
                    alphas = (1-gamma_FW_step)*alphas;
                    alphas(cache_towards_index) = alphas(cache_towards_index) + gamma_FW_step;
                case 'away'
                    alphas = (1+gamma_away_step)*alphas;
                    alphas(cache_away_index) = alphas(cache_away_index) - gamma_away_step;
                case 'pairwise'
                    alphas(cache_away_index) = alphas(cache_away_index) - gamma_pairwise_step;
                    alphas(cache_towards_index) = alphas(cache_towards_index) + gamma_pairwise_step;
                otherwise
                    error('ERROR: internal error: unknown type of step!')
            end
            
            % check if there is a need to renormalize
            need_renormalize_primal = false;
            % sanity check for drop step:
            if is_drop_step
                if ~isequal(best_step_type, 'FW')
                    if alphas(cache_away_index) > 100*eps
                        fprintf('CAUTION: drop step assertion of %s fails on %f, gamma value = %f\n', mfilename, alphas(cache_away_index), gamma);
                        
                        % need to renormalize alphas
                        roundingError = alphas(cache_away_index);
                        alphas(cache_away_index) = 0;
                        renormalizationFactor = (sum(alphas)+roundingError)/sum(alphas);
                        alphas = alphas*renormalizationFactor;
                        need_renormalize_primal = true;
                    else
                        % if rounding is small, can get away with just zeroing out alpha away
                        alphas(cache_away_index) = 0;
                    end
                else
                    if abs(alphas(cache_towards_index)-1) > 100*eps || any(alphas( [1 : cache_towards_index-1, cache_towards_index+1 : end] ) > 100*eps)
                        need_renormalize_primal = true;
                    end
                    alphas = 0 * alphas;
                    alphas(cache_towards_index) = 1;
                end
            end
            if cache{i}.alpha_sum_to_one && abs(sum(alphas)-1)>10*eps
                % need to renormalize alpha, otherwise assertion get failed
                alphas = alphas / sum(alphas);
                need_renormalize_primal = true;
            end
            if need_renormalize_primal
                % after the renormalization of alphas quite heavy renormalization of primal variables is needed
                [ psi_all, loss_all ] = cache_get_entry_by_index( cache{i}, 1 : cache{i}.cache_size );
                model.wMat(:, i) = (1/(lambda*n))*( psi_all*alphas(:) );
                model.ellMat(i) = (1/n) * ( (loss_all(:))'*alphas(:) );
                if ~isempty(positivity)
                    model.v = sum( model.wMat, 2 );
                    % check for sparsity (v should be full)
                    if issparse(model.v)
                        model.v = full(model.v);
                    end
                    model.w   = model.v;
                    model.w(param.positivity==1) = subplus(model.v(param.positivity==1));
                else
                    model.w = sum( model.wMat, 2 );
                    % check for sparsity (w should be full)
                    if issparse(model.w)
                        model.w = full(model.w);
                    end
                end
                model.ell = sum( model.ellMat );
            end
            
            % save the update of the dual variables
            cache{i} = cache_set_dual_vars( cache{i}, alphas );
            
            if options.doCacheConsistencyCheck
                if ~cache_check_consistency_full( cache{i} )
                    fprintf('ERROR in %s: Full cache consistency check did not pass. You definitely should debug this!\n', mfilename);
                end
            else
                if ~cache_check_consistency( cache{i} )
                    fprintf('ERROR in %s: Fast cache consistency check did not pass. You definitely should debug this!\n', mfilename);
                end
            end
        end
        
        if options.use_away_steps || options.use_pairwise_steps
            %% sanity check on alpha variables:
            % we can do this check if active elemets are never removed from cache - this must be the case for pairwise FW and for FW with away steps
            active_I = ~isnan(cache{i}.l_vec);
            if abs(model.ellMat(i) - 1/n*cache{i}.l_vec(active_I)'*cache{i}.alpha_vec(active_I)) >= 1e-12
                fprintf('CAUTION: cache assertion 1 of %s failed on %f\n', mfilename, abs(model.ellMat(i) - 1/n*cache{i}.l_vec(active_I)'*cache{i}.alpha_vec(active_I)) );
            end
            if norm(model.wMat(:,i) - 1/(lambda*n)*double(cache{i}.AFeat(:,active_I))*cache{i}.alpha_vec(active_I)) >= 1e-12
                fprintf('CAUTION: cache assertion 2 of %s failed on %f\n', mfilename, norm(model.wMat(:,i) - 1/(lambda*n)*cache{i}.AFeat(:,active_I)*cache{i}.alpha_vec(active_I)) );
            end
            if abs(sum(cache{i}.alpha_vec(active_I))-1) >= 1e-12
                fprintf('CAUTION: cache assertion 3 of %s failed on %f\n', mfilename, abs(sum(cache{i}.alpha_vec(active_I))-1) );
            end
        end
    end
end

end % onePass_BCFW_hybrid