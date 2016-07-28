function [model, gaps, cache] = onePass_FW(param, options, model, cache)
% [model, gaps, cache] = onePass_FW(param, options, model, cache)
%
% onePass_FW parforms a batch FW step.
% Designed to be called from solver_BCFW_hybrid.m


%% get problem description
phi = param.featureFn; % for \phi(x,y) feature mapping
loss = param.lossFn; % for L(ytruth, ypred) loss function
oracle = param.oracleFn;
n = param.n; % number of training examples
d = param.d; % dimension of feature mapping
using_sparse_features = isfield(param, 'using_sparse_features') && param.using_sparse_features;
lambda = options.lambda;

if ~exist('model', 'var')
    model = initialize_model_zeros(param);
end

% do removeFromActiveSetIfFullCache only if there no away steps and no pairwise steps
removeFromActiveSetIfFullCache = ~options.use_away_steps && ~options.use_pairwise_steps;

%% check the positivity constraints
if isfield(param, 'positivity')
    positivity = true;
    % check that the constraint are respected
    if ~all(model.w(param.positivity==1)>=0)
        fprintf('CAUTION (%s): The model.w does not respect %d positivity constraints\n', mfilename, sum(model.w(param.positivity==1)<0));
    end
    % need to reconstruct the non truncated version of v from the modelInit
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
    positivity = false;
end

%% Initialize models and losses corresponding to the corner atom 
if using_sparse_features
    w_s = sparse(d,1);
    wsMat = sparse(d,n);
else
    w_s = zeros(d,1);
    wsMat = zeros(d,n);
end
ell_s = 0;
ell_s_mat = zeros(n,1);

% if this step is combine with the away/pairwise steps the dual variables have to bee updated even if the cache is not used
if isfield(options, 'update_dual_vars');
    flag_update_dual_vars = options.update_dual_vars;    
else
    flag_update_dual_vars = options.useCache;
end

if flag_update_dual_vars
    % create the yhash table
    yhash = cell(n, 1);
end

%% single pass through data
for i = 1:n
    
    % solve the loss-augmented inference for point i
    ystar_i = oracle(param, model, param.patterns{i}, param.labels{i});
    
    % define the update quantities:
    % [note that lambda*w_s is subgradient of 1/n*H_i(w) ]
    % psi_i(y) := phi(x_i,y_i) - phi(x_i, y)
    psi_i =   phi(param, param.patterns{i}, param.labels{i}) ...
        - phi(param, param.patterns{i}, ystar_i);
    
    wsMat(:,i) = 1/(n*lambda) * psi_i;
    ell_s_mat(i) = loss(param, param.labels{i}, ystar_i)/n;
    
    % sanity check, if this assertion fails, probably there is a bug in the
    % maxOracle or in the featuremap
    % [i.e. the returned loss-augmented inference solution is worse
    % than just using the true label (which has 0 score); this could be
    % caused because a) bug in code; b) the approximate oracle is too
    % bad; c) the true label is *not feasible* ]
    if((n*ell_s_mat(i) - model.w'*psi_i) < -1e-12)
        fprintf('CAUTION: inside %s, probably bug in maxOracle or in the featuremap: %f\n', mfilename, (n*ell_s_mat(i) - model.w'*psi_i));
        ystar_i = param.labels{i};
        psi_i = zeros(size(psi_i));
        wsMat(:,i) = 1/(n*lambda) * psi_i;
        ell_s_mat(i) = loss(param, param.labels{i}, ystar_i)/n;
    end
    
    w_s = w_s + wsMat(:,i);
    ell_s = ell_s + ell_s_mat(i);
    
    if flag_update_dual_vars
        %update the cash values
        yhash{i} = param.hashFn( ystar_i );
    end
end


%% compute duality gap:
gaps = lambda*((model.wMat - wsMat)'*model.w) - (model.ellMat - ell_s_mat);

%% get the step size
if positivity
    gap = lambda*(model.w'*(model.v - w_s)) - (model.ell - ell_s);
    gamma_opt = gap / (lambda*( (model.v - w_s)'*(model.v - w_s) +eps));
else
    gap = lambda*(model.w'*(model.w - w_s)) - (model.ell - ell_s);
    gamma_opt = gap / (lambda*( (model.w - w_s)'*(model.w - w_s) +eps));
end

gamma = max(0, min(1,gamma_opt)); % truncate on [0,1]

if flag_update_dual_vars
    if gamma > 0
        for i = 1 : n
            loss_augmented_scores = cache_get_scores( cache{i}, model.w );
            
            % try to get object from cache, if it is not there it will be added
            [ cache{i}, cache_towards_index ] = cache_get_entry_by_hash( cache{i}, yhash{i}, wsMat(:,i)*(n*lambda),ell_s_mat(i)*n, ...
                loss_augmented_scores, options.maxCacheSize, removeFromActiveSetIfFullCache );
            
            % update the alpha variables:
            alphas = cache_get_dual_vars( cache{i} );
            alphas = (1-gamma)*alphas;
            alphas(cache_towards_index) = alphas(cache_towards_index) + gamma;
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
    end
end

%sanity check
if(abs(sum(gaps) - gap)) > 1E-6
    fprintf('\n\tCaution: Something wrong with %s...\n', mfilename);
end

%% do the step
if gamma > 0
    % finally update the weights and ell variables
    % check for positivity constraints
    if positivity
        model.v = (1-gamma)*model.v   + gamma*w_s;
        % truncate w (optimize over \beta)
        model.w   = model.v;
        model.w(param.positivity==1) = subplus(model.v(param.positivity==1));
    else
        model.w   = (1-gamma)*model.w   + gamma*w_s;
    end
    model.ell = (1-gamma)*model.ell + gamma*ell_s;
    
    %update wMat and ellMat using global gamma
    model.wMat   = (1-gamma)*model.wMat + gamma*wsMat;
    model.ellMat = (1-gamma)*model.ellMat + gamma*ell_s_mat; %check
end

end % onePass_FW
