function [model, gap_vec_heuristic, num_passes, progress] = solver_multiLambda_BCFW_hybrid( param, options )
% [model, gap_vec_heuristic, num_passes, progress] = solver_multiLambda_BCFW_hybrid( param, options )
%
% solver_multiLambda_BCFW_hybrid solves the SSVM problem for multiple values of regularization parameter lambda.
% Supported regimes: grid search with/without warm start and eps-approximate/heuristic regularization path.
%
% The description of the methods is provided in the paper
%
%  [A]  Anton Osokin, Jean-Baptiste Alayrac, Isabella Lukasewitz, Puneet K. Dokania, Simon Lacoste-Julien
%       Minding the Gaps for Block Frank-Wolfe Optimization of Structured SVMs,
%       International Conference on Machine Learning, 2016
%
% Please, cite the paper in any resulting publications.
%
% Fuction solver_multiLambda_BCFW_hybrid relies on solver_BCFW_hybrid.m as an optimization routine for one lambda.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: param, options
%   
%  param: a structure describing the problem. Identical to the on eof solver_BCFW_hybrid.m

%  options:    (an optional) structure with some of the following fields to customize the optimization algorithm:
%
%               Fields sample, gapDegree, gap_check, stepType, useCache, cacheNu, cacheFactor, maxCacheSize, 
%                        doCacheConsistencyCheck, rand_seed, quit_passes_heuristic_gap, quit_passes_heuristic_gap_eps_multiplyer,
%                        logging_level, do_batch_step
%               control the solver for one value of lambda: solver_BCFW_hybrid.m
%
%               regularization_path --controls whether to compute the regularization path or do the grid search
%                                     (default: false) - do grid search
%               
%               true_gap_when_converged -- flag saying whether to demand computation of the true gap at convergence
%                                          (default: true) - eps-approximate path or grid search with guarantees
%
%               check_lambda_change --flag saying whether to do appropriate checks when lambda is changed, significantly slows down the code
%                                     (default: false)
%
%%%%%%%%%% GRID SEARCH
%
%                lambda_values --vector of lambdas defining the grid
%                                (default: 10.^(4: -1: -3)
%
%                gap_threshold --accuracy, until which to optimize for each lambda
%                                (default: 0.1)
%
%                warm_start_type -- type of warm start to use: 'keep_primal' or 'keep_dual' or 'none';
%                                   (default: 'keep_primal')
% 
%%%%%%%%%% REGULARIZATION PATH
%
%                regularization_path_min_lambda -- minimum value of lambda determining the last breakpoint
%                                                  The method will quit earlier if the stopping criterion is hit.
%                                                  (default: 1e-5)
%
%                regularization_path_eps --eps parameter for eps-approximate path 
%                                          (default: 1e-1)
%
%                regularization_path_a   --kappa parameter from the paper. The internal SSVM server is
%                                          run with kappa*eps gap stopping criterion
%                                          (default: 0.9)
%
%                regularization_path_b   --when doing the induction step to get the new breakpoint we can be 
%                                          less conservative of more conservative. 
%                                          (default: 1 - options.regularization_path_a)
%
%%%%%%%%%% COMPUTATIONAL BUDGET (whole process)
%
%                num_passes   -- max number of passes through data
%                                (default: 200000) - effectively unlimites
%                time_budget  -- max running time
%                                (default: 60*24) - 24 hours
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: model, gap_vec_heuristic, num_passes, progress, cache, exact_gap_flag
%
%   model     --model.w contains the obtained parameters w (column vector)
%               model.ell contains b'*alpha which is useful to compute duality gap, etc.
%               model.wMat contains one parameter vector w_i per training object
%               model.ellMat contains one value \ell_i per training object
%               model.v is present only when there are come positivity constraints (untruncated version of model.w)
%
%   gap_vec_heuristic  -- the estimates of the block gaps obtained at the end of the method.
%                         If the method has converged and options.true_gap_when_converged == true
%                         the estimates equal the exact block gaps.
%
%   num_passes -- number of effective passes over the dataset performed by the method
%
%   progress  --logged information about the run of the method. The ammount of information are determined by options.logging_level
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Anton Osokin, Jean-Baptiste Alayrac, Simon Lacoste-Julien
% Project web page: http://www.di.ens.fr/sierra/research/gapBCFW
% Code: https://github.com/aosokin/gapBCFW



%% parse the options
options_default = defaultOptions;
if exist('options', 'var')
    options = processOptions(options, options_default);
else
    options = options_default;
end
fprintf('Running %s on %d examples.\n', mfilename, numel(param.patterns));
switch options.stepType
    case 0
        options.use_FW_steps = true;
        options.use_away_steps = true;
        options.use_pairwise_steps = true;
        fprintf('Launching the combined block coordinate method: FW, away and pairwise steps\n');
    case 1
        options.use_FW_steps = true;
        options.use_away_steps = false;
        options.use_pairwise_steps = false;
        fprintf('Launching block coordinate FW (BCFW)\n');
    case 2
        options.use_FW_steps = false;
        options.use_away_steps = false;
        options.use_pairwise_steps = true;
        fprintf('Launching block coordinate pairwise FW (BC-P-FW)\n');
    case 3
        options.use_FW_steps = true;
        options.use_away_steps = true;
        options.use_pairwise_steps = false;
        fprintf('Launching block coordinate FW with away steps (BC-A-FW)\n');
    otherwise
        error([mfilename, ': unknown execution mode'])
end
% if cache is used or we are going to do pairwise or away steps we will need to create and maintain dual variables
options.update_dual_vars = options.useCache || options.use_away_steps || options.use_pairwise_steps;
% values of lambda need to be decreasing
options.lambda_values = sort(options.lambda_values, 'descend');
% the following flag must always be false when steps with away corners are used
% otherwise solver_multiLambda_BCFW_hybrid.m will be outputting dual variables inconsistent with the primal
options.output_model_before_batch_step = false;
% when we do regularization path we optimize at each lambda up to A*eps gap
if options.regularization_path
    options.gap_threshold = options.regularization_path_a * options.regularization_path_eps;
end
% the regularization path mode
% when options.regularization_path == true parameters options.lambda_values and options.gap_threshold are ignored
% type of warm start; has to be 'keep_primal' for regularization path
if options.regularization_path
    options.warm_start_type = 'keep_primal'; 
end

fprintf('The options are as follows:\n');
options


%% start the timer
tStart = tic();

%% fix the seed
init_rand_seed = options.rand_seed;
rng('default');
rng(init_rand_seed);

%% check that the loss on ground truth itself is zero:
% param.lossFn(param, label_gt, label_gt) == 0
% otherwise some of the equations below are incorrect
for i_object = 1 : param.n
    label_gt = param.labels{i_object};
    loss_gt = param.lossFn(param, label_gt, label_gt);
    if loss_gt>1e-12
        error(['On object ', num2str(i_object), ' loss netween GT and GT is non zero. This case is currently not supported.']);
    end
end

%% initialization
if ~options.regularization_path
    % without the reg. path we do regular initialization with zeros
    if options.update_dual_vars
        % create cache for zero-initialized model
        [cache, model] = cache_initilize_model( param, options.maxCacheSize );
    else
        cache = [];
        if ~exist('model', 'var')
            model = initialize_model_zeros(param);
        end
    end
    change_lambda_flag = false;
    lambda_previous = inf;
    exact_gap_flag = false;
    if ~exist('gap_vec_heuristic', 'var')
        gap_vec_heuristic = inf(param.n,1);
    end
else
    % in case of the regularization path we need to find lambda, primal and dual variables, such that gap is small enough
    fprintf('Initializing the regularizarion path\n');
    [options.lambda_values, model, gap_vec_heuristic, cache] = initialize_regularization_path( param, options.gap_threshold, options.update_dual_vars, options.maxCacheSize );
    change_lambda_flag = true;
    lambda_previous = options.lambda_values;
    exact_gap_flag = false;
    options.init_gap_value_for_cache = sum(gap_vec_heuristic);
end

%% do the main loop
num_lambdas = numel(options.lambda_values);
progress = cell(0, 1);
i_lambda = 0;
time_budget_full = options.time_budget;
num_passes = 0;
while (~options.regularization_path && i_lambda < numel(options.lambda_values)) || ...
        (options.regularization_path && (i_lambda==0 || options.lambda_values( i_lambda ) > options.regularization_path_min_lambda))
    %% init pass for the new value of lambda
    i_lambda = i_lambda + 1;
    options.rand_seed = init_rand_seed + i_lambda - 1;
    tStartLambda = tic();
    
    %% get the new value of lambda
    % compute the vector update for the gaps, it is used for the reg. path and later to update the gap estimate
    gap_vec_heuristic_update = model.ellMat - lambda_previous*(model.wMat'*model.w(:));
        
    if ~options.regularization_path
        options.lambda = options.lambda_values( i_lambda );
        fprintf('Processing lambda %f: %d of %d\n', options.lambda, i_lambda, num_lambdas);
    else
        % current gap estimate is supposed to be smaller than options.regularization_path_a * options.regularization_path_eps
        gap_heuristic = sum( gap_vec_heuristic );
        if gap_heuristic > options.regularization_path_a * options.regularization_path_eps + 1e-12
            fprintf('CAUTION (%s): current gap estimate=%f is greater than A*eps= %f, something probably went wrong.\n', mfilename, gap_heuristic, options.regularization_path_a*options.regularization_path_eps);
        end
        gap_heuristic_update = sum( gap_vec_heuristic_update );

        % we can do larger steps taking into account the current estimates of the gaps
        max_allowed_gap_update = (options.regularization_path_b + options.regularization_path_a) * options.regularization_path_eps - gap_heuristic;
        
        % small update stopping criterion
        if gap_heuristic_update <= max_allowed_gap_update
            % this case happens when gap_heuristic_update <= B * eps
            % it means that for any lambdas smaller than the current one new gap will be smaller than (A+B)*eps
            fprintf('Gap update is smaller than B*eps. Terminating.\n');
            break;
        end

        % computing the maximum change of lambda that we can do
        lambda_ratio = 1 - max_allowed_gap_update / gap_heuristic_update;
        if lambda_ratio > 1
            error('CAUTION: next lambda is going to be larger than the previous one. Heuristic gap estimate is negative! Something is wrong!');
        end
        if lambda_ratio <= 0
            error('CAUTION: next lambda is going to be negative. Something is wrong!');
        end
        
        options.lambda = lambda_ratio * lambda_previous;
        options.lambda_values = [options.lambda_values; options.lambda];
        fprintf('Processing lambda %f: %d of the regularization path, min lambda: %f\n', options.lambda, i_lambda, options.regularization_path_min_lambda);
    end
    
    %% update model and cache when changing lambda
    can_do_warm_start = true;
    if isequal( options.warm_start_type, 'keep_primal' )
    if change_lambda_flag
        % update the model when changing lambda
        lambda_ratio = options.lambda / lambda_previous;
        % to keep result of the oracle the same, we keep primal variables model.w, model.wMat
        
        % update the gap values
        % if the old gap value was exact than the new gap value will be exact as well    
        if ~exact_gap_flag
            % if the old gap estimates were not correct we do not decrease the gap to be conservative
            gap_vec_heuristic_update( gap_vec_heuristic_update < 0 ) = 0;
        end
        gap_vec_heuristic = gap_vec_heuristic + (1-lambda_ratio)*gap_vec_heuristic_update;
        
        % update values model.ell, model.ellMat
        model.ell = model.ell * lambda_ratio;
        model.ellMat = model.ellMat * lambda_ratio;
        
        % the dual variables need to be changed if maintained
        if options.update_dual_vars
            for i_object = 1 : param.n
                % find index of the GT label
                label_gt = param.labels{i_object};
                yhash_i_gt = param.hashFn( label_gt );
                [cache{i_object}, cache_gt_index ] = cache_get_entry_by_hash( cache{i_object}, yhash_i_gt, zeros(param.d, 1), param.lossFn(param, label_gt, label_gt), ...
                    [], options.maxCacheSize, false );
                alphas = cache_get_dual_vars( cache{i_object} );
                
                if ~cache{i_object}.alpha_sum_to_one
                    fprintf('CAUTION (%s): Alphas of object %d in the cache do not sum up to one: can not reinitialize for another lambda.', mfilename, i_object);
                    can_do_warm_start = false;
                end
                
                % all alphas corresponding to non-groundtruth labelings are multiplies by ratios lambda_new / lambda_old
                % the rest of the mass goes into the ground-truth lambda
                alphas = alphas * lambda_ratio;
                alphas(cache_gt_index) = 0;
                alphas(cache_gt_index) = 1 - sum(alphas);
                
                cache{i_object} = cache_set_dual_vars( cache{i_object}, alphas );
            end
        end
    end
    
    %% check model and cache after the update
    if options.check_lambda_change
        if options.update_dual_vars
            % check if the model reconstructed from the cache is eqauivalent to the maintained one
            if options.do_batch_step && options.output_model_before_batch_step
                % when options.do_batch_step==true and options.output_model_before_batch_step==true
                % checking the gap does not make sense, because it is outdated after the batch FW step
                % only methods not relying on the dual variables should be used in this regime
                fprintf('Skipping the check of the model in %s because it is outdated in this regime.\n', mfilename);
                if options.use_away_steps || options.use_pairwise_steps
                    fprintf('CAUTION (%s): away and pairwise steps will be incorrect in this regime because the dual variables are not consistent with the primal ones. Use at your own risk!\n', mfilename);
                end
            else
                [ cache, model_dual ] = cache_initilize_model( param, options.maxCacheSize, cache, [], options.lambda );
                if ~isequal_models(model, model_dual)
                    fprintf('CAUTION (%s): the update to lambda %f leads to inconsistent model\n', mfilename, options.lambda);
                end
            end
        end
        if exact_gap_flag
            % checking the gap values makes sense only the maintained gaps equal the true gaps
            [~, gap_vec_check] = duality_gap_vec( param, param.oracleFn, model, options.lambda, model.wMat, model.ellMat );
            if any( abs(gap_vec_heuristic(:) - gap_vec_check(:)) > 1e-8 )
                fprintf('CAUTION (%s): gaps reconstructed for lambda %f do not equal the maintained ones, max diff: %f\n', mfilename, options.lambda, max(abs(gap_vec_heuristic(:) - gap_vec_check(:) )) );
            end
            
        end
    end
    elseif isequal( options.warm_start_type, 'keep_dual' )
        if ~options.update_dual_vars
            error('Cannot warm start from dual variables, because dual variables are not initialized.');
        end
        % construct primal model from dual variables stored in cache
        [ cache, model ] = cache_initilize_model( param, options.maxCacheSize, cache, [], options.lambda );
        % with this type of warm start we have not control of what is happening with the gaps
        % we can either set the gaps to infinity or do a batch pass to update them
        
        % [~, gap_vec_heuristic] = duality_gap_vec( param, param.oracleFn, model, options.lambda, model.wMat, model.ellMat );
        % exact_gap_flag = true;
        
        gap_vec_heuristic = inf(param.n,1);
        exact_gap_flag = false;
    elseif isequal( options.warm_start_type, 'none' )
        can_do_warm_start = false;
    else
        error('Unknown warm start type!');
    end
    
    %% run the solver on the new value of lambda
    % update the time budget for the solver:
    options.time_budget = time_budget_full - toc( tStart )/60;
    if  exact_gap_flag
        options.init_gap_value_for_cache = sum(gap_vec_heuristic);
    else
        % in this case we just turn off global criterion for cache
        options.init_gap_value_for_cache = 0;
    end
    if can_do_warm_start
        [ model, gap_vec_heuristic, curNumPasses, progress{i_lambda}, cache, exact_gap_flag ] = solver_BCFW_hybrid(param, options, model, gap_vec_heuristic, cache);
    else
        [ model, gap_vec_heuristic, curNumPasses, progress{i_lambda}, cache, exact_gap_flag ] = solver_BCFW_hybrid(param, options);
    end
    num_passes = num_passes + curNumPasses;
    
    %% time budget stopping criterion
    if toc(tStart)/60 > time_budget_full
        % time budget fully used
        fprintf('Spent %fs on lambda %f; did not converge\n', toc(tStartLambda), options.lambda);
        fprintf('Time budget was fully used. Outputting the part that was computed\n');
        break;
    end
    
    %% final touches
    lambda_previous = options.lambda;
    change_lambda_flag = true;
    fprintf('Spent %fs on lambda %f\n', toc(tStartLambda), options.lambda);
end

%% final steps
fprintf('Spent %fs for %d values of lambda\n', toc(tStart), num_lambdas )
end % solverBCFWH

function equal = isequal_models(model1, model2, accuracy)
if ~exist('accuracy', 'var') || isempty(accuracy)
    accuracy = 1e-12;
end
equal = true;
if any( abs(model1.w(:) - model2.w(:)) > accuracy )
    fprintf('CAUTION (%s): field <w> of the two models has max difference of %f\n', mfilename, max(abs(model1.w(:) - model2.w(:))) );
    equal= false;
end
if any( abs(model1.ell(:) - model2.ell(:)) > accuracy )
    fprintf('CAUTION (%s): field <ell> of the two models has max difference of %f\n', mfilename, max(abs(model1.ell(:) - model2.ell(:))) );
    equal= false;
end
if any( abs(model1.wMat(:) - model2.wMat(:)) > accuracy )
    fprintf('CAUTION (%s): field <wMat> of the two models has max difference of %f\n', mfilename, max(abs(model1.wMat(:) - model2.wMat(:))) );
    equal= false;
end
if any( abs(model1.ellMat(:) - model2.ellMat(:)) > accuracy )
    fprintf('CAUTION (%s): field <ellMat> of the two models has max difference of %f\n', mfilename, max(abs(model1.ellMat(:) - model2.ellMat(:))) );
    equal= false;
end
if isfield(model1, 'v') ~= isfield(model2, 'v')
    fprintf('CAUTION (%s): the two models are inconsistent in terms of the support of positivity constraints\n', mfilename);
    equal= false;
end
if isfield(model1, 'v') && isfield(model2, 'v') && any( abs(model1.v(:) - model2.v(:)) > accuracy )
    fprintf('CAUTION (%s): field <v> of the two models has max difference of %f\n', mfilename, max(abs(model1.v(:) - model2.v(:))) );
    equal= false;
end
end

function options = defaultOptions
options = struct;

% for grid search
options.lambda_values = 10.^(4: -1: -3);
options.gap_threshold = 0.1;

% stopping parameters. joint budget for whole solver_multiLambda_BCFW_hybrid.m
options.num_passes = 200000; % max number of passes through data
options.time_budget = 60*24;

% flag saying whether to demand computation of the true gap at the convergence
options.true_gap_when_converged = true;

% flag saying whether to do appropriate checks when lambda is changed, significantly slows down the code
options.check_lambda_change = false;

% type of warm start; has to be 'keep_primal' for regularization path
options.warm_start_type = 'keep_primal'; % 'keep_primal' or 'keep_dual' or 'none';

% the regularization path mode
% when options.regularization_path == true parameters options.lambda_values and options.gap_threshold are ignored
options.regularization_path = false;
options.regularization_path_eps = 1e-1;
options.regularization_path_a = 0.9;
options.regularization_path_b = 1 - options.regularization_path_a;
options.regularization_path_min_lambda = 1e-5;

%% parameters from solver_BCFW_hybrid.m
% gap sampling: 'uniform' or 'perm' or 'gap' or 'maxGap'
options.sample = 'gap';
options.gapDegree = 1.0;

% cache options
options.useCache = true;
options.cacheNu = 0.01;
options.cacheFactor = 0.25;
options.maxCacheSize = 100;
options.doCacheConsistencyCheck = false; % simple check is always done anyway; if true slow exhaustive check is done

% choose in which mode to run the combined solver
% 0 - all steps
% 1 - only FW
% 2 - pairwise
% 3 - FW and away
options.stepType = 2;

% do batch FW step or just compute the duality gap
options.do_batch_step = false;

% other parameters
options.gap_check = 10; % how often to compute the true gap
options.rand_seed = 1; % random seed

% flag to exit onePass_BCFW_hybrid.m based on the heuristic gaps computed inside
% this flag affects the method, but convergence guarantees at the end at still exact
options.quit_passes_heuristic_gap = true;
options.quit_passes_heuristic_gap_eps_multiplyer = 0.8; % quit BC passes if heuristic gap is this factor smaller than gap_threshold;  < 1 for underetimated gaps; > 1 for overestimated gaps

% the logging level, how many things to store?
% - 0: no logging at all
% - 1: save models after each pass, and some data statistics, should not cause huge memory and speed overheads
% - 2: save everything reasonable, including O(n) numbers per each gap check pass
% - 3: crazy logging, save everything possible per point at each pass; huge computational overheads, use this only for specific plots
options.logging_level = 1;

end % defaultOptions


