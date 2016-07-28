function [model, gap_vec_heuristic, num_passes, progress, cache, exact_gap_flag] = solver_BCFW_hybrid(param, options, model, gap_vec_heuristic, cache)
% [model, gap_vec_heuristic, num_passes, progress, cache, exact_gap_flag] = solver_BCFW_hybrid(param, options, model, gap_vec_heuristic, cache)
%
% Solves the structured support vector machine (SSVM) using different variants of the block-coordinate
% Frank-Wolfe (BCFW) algorithm. See
%
% [A]   Anton Osokin, Jean-Baptiste Alayrac, Isabella Lukasewitz, Puneet K. Dokania, Simon Lacoste-Julien
%       Minding the Gaps for Block Frank-Wolfe Optimization of Structured SVMs,
%       International Conference on Machine Learning, 2016
%
% [B]   Simon Lacoste-Julien, Martin Jaggi, Mark Schmidt, Patrick Pletscher
%       Block-Coordinate Frank-Wolfe Optimization for Structural SVMs,
%       International Conference on Machine Learning, 2013
%
% for more details. If you use this code please cite these papers in any resulting publications.
%
% This code was originally forked from BCFWstruct (https://github.com/ppletscher/BCFWstruct) 
% by Simon Lacoste-Julien, Martin Jaggi, Mark Schmidt, Patrick Pletscher.
% The interface is similar.
%
%
% The SSVM has the form
% min_{w} 0.5*\lambda*||w||^2+ 1/n*\sum_{i=1}^n H_i(w) [ see (1) in paper [A] ]
%   where H_i(w) is the structured hinge loss on the i^th example:
%         H_i(w) = max_{y in Y} L(y_i,y) - <w, \psi_i(y)> [ (2) in paper [A] ]
%
% We use a calling interface very similar to version 1.1 of svm-struct-matlab developped by Andrea Vedaldi 
% (see vedaldi/code/svm-struct-matlab.html). svm-struct-matlab is a Matlab wrapper interface to the widely used
% SVM^struct code by Thorsten Joachims (http://svmlight.joachims.org/svm_struct.html) which implements a cutting plane algorithm.
% See paper [B] for experiments showing that the convergence per oracle call is much slower for SVM^struct than BCFW.
%
% If your code was using
%    model = svm_struct_learn(command_args,param)
% or
%    model = solverBCFW(param, options)
% you can replace it with
%    model = solver_BCFW_hybrid(param, options)
%
% All the other inputs/outputs are optional.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: param, options, model, gap_vec_heuristic, cache
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  param: a structure describing the problem. The struture can contain the dataset and is never copied. 
%            Extra fields can be added to param, e.g., to pass some extra data into the oracle and feature map.
%            param contains the following fields:      
%
%     n         -- number of training objects
%
%     d         -- number of trainable parameters, dimensionality of the feature map 
%
%     patterns  -- patterns (x_i)
%         A cell array of patterns (x_i). The entries can have any nature (they can just be indexes 
%         of the actual data, for example).
%
%     labels    -- labels (y_i)
%         A cell array of labels (y_i). The entries can have any nature.
%
%     lossFn    -- loss function callback
%         A handle to the loss function L(ytruth, ypredict) defined for your problem. This function 
%         should have a signature of the form:
%             scalar_output = loss(param, ytruth, ypredict)
%         It will be given an input ytruth, a ground truth label; ypredict, a prediction label; 
%         and param, the same structure passed to solver_BCFW_hybrid.
%
%     oracleFn  -- loss-augmented decoding callback
%         A handle to the 'maximization oracle' (equation (2) in paper [A]) which solves the loss-augmented 
%         decoding problem. This function should have a signature of the form:
%             ypredict = decode(param, model, x, y)
%         where x is an input pattern, y is its ground truth label, param is the input param structure 
%         to solver_BCFW_hybrid and model is the current model structure (the main field is model.w which contains
%         the parameter vector).
%
%     featureFn -- feature map callback
%         A handle to the feature map function \phi(x,y). This function should have a signature of the form:
%           phi_vector = feature(param, x, y) --- COLUMN VECTOR
%         where x is an input pattern, y is an input label, and param is the usual input param structure. 
%         The output should be a vector of *fixed* dimension d which is the same across all calls to the function. 
%         The parameter vector w will have the same dimension as this feature map. In our current
%         implementation, w is sparse if phi_vector is sparse.
%
%     hashFn   -- encodes the output with a unique string, required when using caching, away or pairwise steps.
%         Takes a label as input and outputs a string of characters.
%         This is used to create the dual variables and identity them uniquely.
%
%     positivity  -- vector of length d, with ones indicating that the corresponding parameter always has to be positive
%           (optional, default: zeros(d,1) - no constraints)
%
%     using_sparse_features -- flag saying whether feature operations need to be done with sparse vectors
%           (optional, default: false)
%
%     cache_type_single -- flag saying whether the entries of cache need to be stored in single instead of double
%           (optional, default: false)
%           incompatible with using_sparse_features == true
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  options:    (an optional) structure with some of the following fields to customize the optimization algorithm:
%
%%%%%%%%%% REGULARIZATION
%
%   lambda    --The regularization constant (default: 1/n).
%
%%%%%%%%%% GAP SAMPLING, See paper [A], Section 3.1 for details
%
%   sample    --Sampling strategy for example index. Options:
%                   'uniform' - sampling indices from uniform distribution on {1, ..., n}
%                   'perm' - sampling indices according to a random permutation ant each pass
%                           [Note that convergence-rate proofs of [B] only hold for uniform
%                           sampling, not for a random permutation.]
%                   'gap' - sample according to gap, i.e., proportional to g_i^p, where g_i are the block estimates
%                           p is specified by optional parameter options.gapDegree (default: 1)
%                   'maxGap' - deterministic gap sampling, picking the object with maximum gap estimate.
%                              CAUTION: leads to severe staleness effect.
%               (default: 'gap')
%
%   gap_check --Gives the number of passes through data between each check of the duality gap (put inf to turn it off).
%               Every gap_check passes through the data, the algorithm does a batch pass through the data to compute
%               the duality gap or to do a batch FW step. This slows down the code roughly by a factor 
%               of ( 1 + 1/(gap_check+1) ).
%               This is the parameter controlling exploitation-staleness trade-off
%               (default: 10)
%
%%%%%%%%%% TYPES OF FW STEPS, See paper [A], Section 3.2 for details
%
%   stepType  --choose steps to use in the solver
%                   0 - all steps
%                   1 - only regular FW
%                   2 - pairwise 
%                   3 - FW and away
%               (default: 2)
%
%%%%%%%%%% CACHING, See paper [A], Section 3.3 for details
%
%   useCache  --flag saying whether to use cache or not.
%               Criterion for the cache hit (not calling the oracle): 
%                  gap_cache >= max(g_i*options.cacheFactor, options.cacheNu/n*fullGapEstimate)
%               (default: true)
%
%   cacheNu   --parameter for global cache hit criterion.
%               Larger values lead to more oracle calls, but have better constant in the convergence guarantee.
%               (default: 0.01)
%
%   cacheFactor   --parameter for local cache hit criterion.
%                   (default: 0.25)
%
%   maxCacheSize  --maximum allowed size of the cache.
%                   (default: 100)
%
%   doCacheConsistencyCheck  --flag controling the consistency checks of the cache.
%                              If false only a very cheap consistency check is done; 
%                              If true slow exhaustive check is done
%                              (default: false) 
%
%%%%%%%%%% STOPPING CRITERIA
%
%   gap_threshold 
%               --Stop the algorithm once the duality gap falls below this threshold. Use gap_check to 
%                 control how often the gap is checked. Note that the default of 0.1 is equivalent
%                 to the criterion used in svm_struct_learn Matlab wrapper.
%                 (default: 0.1)
%
%   num_passes  --Maximum number of passes through the data before the algorithm stops.
%                 Note, this is not the number of EFFECTIVE passes, i.e. cache hits are counted.
%                 (default: 200)
%
%   time_budget --Number of minutes after which the algorithm should terminate.
%                 Useful if the solver is run on a cluster with some runtime limits. 
%                 (default: inf)
%
%%%%%%%%%% DIFFERENT
%
%   rand_seed --Optional seed value for the random number generator.
%               (default: 1)
%
%   do_batch_step --flag, saying whether to do batch FW step or just compute the duality gap
%                   (default: false) - just compute the gap
%
%   true_gap_when_converged --flag controlling whether to demand exact gap guarantee or exit with heuristic gaps
%                             (default: true) - exact gap guarantee   
%
%   logging_level    --the logging level, how many things to store?
%                         - 0: no logging at all
%                         - 0.5: save models after each gap checking pass, should cause minimal memory and speed overheads 
%                         - 1: save models after each pass, and some data statistics, should not cause large memory and speed overheads
%                         - 2: save everything reasonable, including O(n) numbers per each gap check pass (small computational overhead, but significant storage)
%                         - 3: crazy logging, save everything possible per point at each pass; huge computational overheads, use this only for specific plots
%                      (default: 1)
%
%   quit_passes_heuristic_gap  --flag, saying if exit onePass_BCFW_hybrid.m based on the heuristic gaps 
%                                computed inside. This flag affects the method, but convergence guarantees 
%                                at the end are still exact
%                                (default: true) - exit the inner loop heuristically
%
%   quit_passes_heuristic_gap_eps_multiplyer   --parameter, controlling when to exit the loop according to heuristic gaps.
%                                                Less than one to be more conservative.
%                                                (default: 0.8)
%
%   output_model_before_batch_step   --If you actually want to use the output with batch FW steps (e.g. measure test error) 
%                                      we recommend to set this flag to true.  
%                                      We observed that after the batch FW step the primal objective significantly goes up (the gap is still going down).
%                                      However, if you want to use the output dual variables (cache) the following flag must always be false, 
%                                      otherwise solver_BCFW_hybrid.m will be outputting dual variables inconsistent with the primal.
%                                      Does not have any effect when options.do_batch_step == false.
%                                      (Default: false)
%
%   init_gap_value_for_cache   --initialization value of global gap used for the global cache criterion
%                                   inf - always fail
%                                   0 - always pass
%                                Is active until the first batch gap computation pass.   
%                                (default: 1) - is correct if the loss is bounded by one, otherwise heuristics until the first batch pass 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Inputs model, gap_vec_heuristic, and cache are used to warm start the solver.
%       If you do not want to do warm start just leave them empty.
%       The structure of these Inputs matches the corresponding Outputs.
%       The details are provided if the section for Outputs.
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
%   cache --the cache data structure obtained at convergence. Required to do the warm start.
%           cell array, one for each example, containing the structure representing the dual variables:
%                   cached_labels: this is a dictionary (containers.Map) mapping labeling to the index in the cache
%                                  each entry of the dictionary is: (y_key; index_for_feature); 
%                                  y_key is computed using param.hashFn(ylabel)
%                   AFeat: this is a size of cache x d matrix, each row storing the feature vector for a specific output
%                   alpha_vec: each row is entry for alpha variable; 0 when inactive...
%                   l_vec: each row contains the loss for this output
%                   keys: each row gives the key of the corresponding output in the active set...
%                   cache_size: size of the current cache
%                   cache_entries: cache positions, occupied by some elements
%                   alpha_sum_to_one: flag showing whether model thinks that alphas is cache sum to one, i.e. full active set is in the cache
%    
%   exact_gap_flag --a flag showing whether the gap computed at the end of the process is exact
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Anton Osokin, Jean-Baptiste Alayrac, Simon Lacoste-Julien
% Project web page: http://www.di.ens.fr/sierra/research/gapBCFW
% Code: https://github.com/aosokin/gapBCFW


%% parse the options
options_default = defaultOptions(param.n);
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

fprintf('The options are as follows:\n');
options


%% fix the seed
rng('default');
rng(options.rand_seed);

%% init progress report
progress = struct;
if options.logging_level >= 1
    progress.models = cell(0, 0); % saving all the intermediate models
    progress.time = []; % wall-clock time
    progress.numPasses = []; % effective number of passes: number of oracle calls / n
    progress.verifiedGap = []; % exact values of the gap when it is computed (no overhead)
    progress.verifiedGapModelId = []; % index of the model, for which the exact gap is computed
    progress.verifiedGapNumPasses = []; % number of effective passes when the exact gap was computed
    progress.lambda = options.lambda; % save lambda, just in case
    progress.data_pass_info =  cell(0, 0);
    
    if options.logging_level >= 2
        progress.trueGapsPerPoint = cell(0, 0);
        progress.heuristicGapsPerPoint = cell(0, 0);
        if options.useCache
            progress.cacheSizePerPoint =  cell(0, 0);
            progress.activeSetSizePerPoint =  cell(0, 0);
        end
    end
end

tStart = tic();

%% initialize cache
if options.update_dual_vars
    if ~exist('cache', 'var') || isempty( cache )
        % create cache for zero-initialized model
        [cache, model] = cache_initilize_model( param, options.maxCacheSize);
    elseif exist('model', 'var') && ~isempty(model)
        [cache, model] = cache_initilize_model( param, options.maxCacheSize, cache, model );
    else
        [cache, model] = cache_initilize_model( param, options.maxCacheSize, cache, [], options.lambda );
    end
else
    cache = [];
    if ~exist('model', 'var')
        model = initialize_model_zeros(param);
    end
end
if ~exist('gap_vec_heuristic', 'var')
    gap_vec_heuristic = Inf(param.n,1);
end

%% main loop
num_passes = 0;
gap_value_for_cache = options.init_gap_value_for_cache; 
i_pass = 0;
converged_heuristic = false;
exact_gap_flag = false;
while (i_pass < options.num_passes)
    % run bcfw maximum upto gap_check passes
    num_bcfw_pass = 0;
    for i = 1:options.gap_check
        i_pass = i_pass+1;
        fprintf('Pass %d of %d, heuristic gap: %f\n', i_pass, options.num_passes, sum(gap_vec_heuristic));
        [model, gap_vec_heuristic, frac_pass, cache, data_pass_info] = onePass_BCFW_hybrid(param, options, model, gap_vec_heuristic, cache, gap_value_for_cache);
        exact_gap_flag = false;
        num_bcfw_pass = num_bcfw_pass + frac_pass;
        
        if options.logging_level >= 1
            % logging current status
            model_for_saving = struct;
            model_for_saving.w = model.w;
            model_for_saving.ell = model.ell;

            progress.models = [progress.models; {model_for_saving} ];
            progress.time = [progress.time; toc(tStart)];
            progress.numPasses = [progress.numPasses; num_passes+num_bcfw_pass];
            progress.data_pass_info = [progress.data_pass_info; {data_pass_info}];

            if options.logging_level >= 3
                % at each pass save something O(n) + compute the exact values of the gaps
                % for gap samling plots (normally, you really do not want to do this):
                progress.heuristicGapsPerPoint = [progress.heuristicGapsPerPoint; {gap_vec_heuristic}];
                [~, gap_vec_true_recompute] = duality_gap_vec(param, param.oracleFn, model, options.lambda, model.wMat, model.ellMat);
                progress.trueGapsPerPoint = [ progress.trueGapsPerPoint; {gap_vec_true_recompute} ];
                % for caching plots:
                if options.useCache
                    cacheSizePerPoint = zeros(param.n, 1);
                    activeSetSizePerPoint = zeros(param.n, 1);
                    for i_object = 1:param.n
                        cacheSizePerPoint(i_object) = cache{i_object}.cache_size;
                        activeSetSizePerPoint(i_object) = sum( cache{i_object}.alpha_vec(cache{i_object}.cache_entries) > 0 );
                    end
                    progress.cacheSizePerPoint = [progress.cacheSizePerPoint; {cacheSizePerPoint}];
                    progress.activeSetSizePerPoint = [progress.activeSetSizePerPoint; {activeSetSizePerPoint}];
                end
            end
            % end of logging current status
        end
        
        if options.quit_passes_heuristic_gap && (sum(gap_vec_heuristic) <= options.gap_threshold*options.quit_passes_heuristic_gap_eps_multiplyer) % heuristic gap estimates are below the threshold
            converged_heuristic = true;
            break;
        end
        % time-budget exceeded?
        t_elapsed = toc(tStart);
        if (t_elapsed/60 > options.time_budget)
            fprintf('Time budget exceeded.\n');
            break;
        end
        if i_pass > options.num_passes
            fprintf('Maximum number of passes exceeded.\n');
            break;
        end
    end
    
    cache_mem_info = whos('cache');
    all_vars_mem_info = whos;
    progress_mem_info = whos('progress');
    fprintf('Cache size (GB): %f, progress size (GB): %f, all other vars (GB): %f\n', ...
        cache_mem_info.bytes / 1024^3, progress_mem_info.bytes / 1024^3, (sum(cat(1, all_vars_mem_info.bytes))-cache_mem_info.bytes-progress_mem_info.bytes) / 1024^3 );
    
    model_before_batch_step = model;
        
    if converged_heuristic && ~options.true_gap_when_converged
        % in this regime we do not check the true gap after heuristic convergence
        num_passes = num_passes+num_bcfw_pass;
        fprintf('Heuristic duality gap below threshold at pass: %.3f -- stopping!\n', num_passes)
        fprintf('Heuristic gap: %g, gap_threshold: %g\n', sum(gap_vec_heuristic), options.gap_threshold)
        break;
    else
        % normal regime: run duality gap check or batch FW pass
        % if options.do_batch_step == true do the batch FW step; otherwise just compute the duality gap
        if options.do_batch_step
            % run FW pass
            [model, gap_vec_true, cache] = onePass_FW(param, options, model, cache);
            exact_gap_flag = false;
        else
            [~, gap_vec_true] = duality_gap_vec(param, param.oracleFn, model, options.lambda, model.wMat, model.ellMat);
            exact_gap_flag = true;
        end
        gap_true = sum(gap_vec_true);
        num_passes = num_passes+1+num_bcfw_pass;
        i_pass = i_pass+1;
        
        gap_value_for_cache = gap_true; % used for caching
        
        if options.logging_level > 0
            % logging current status
            % logging happening at every pass:
            model_for_saving = struct;
            model_for_saving.w = model.w;
            model_for_saving.ell = model.ell;
            progress.models = [progress.models; {model_for_saving} ];
            progress.time = [progress.time; toc(tStart)];
            progress.numPasses = [progress.numPasses; num_passes];
            data_pass_info = struct;
            data_pass_info.noCacheHit = param.n;
            data_pass_info.globalGapHit = 0;
            data_pass_info.localGapHit = 0;
            data_pass_info.doubleGapHit = 0;
            data_pass_info.numCacheHits = 0;
            data_pass_info.numFwSteps = param.n;
            data_pass_info.numAwaySteps = 0;
            data_pass_info.numPairwiseSteps = 0;
            progress.data_pass_info = [progress.data_pass_info; {data_pass_info}];

            % logging happening at each full gap check:
            progress.verifiedGap = [progress.verifiedGap; gap_true];
            progress.verifiedGapModelId = [progress.verifiedGapModelId; numel(progress.models)-1]; % the verified duality gap check was computed for the previous model, before the FW batch pass
            progress.verifiedGapNumPasses = [progress.verifiedGapNumPasses; num_passes];
            
            if options.logging_level >= 2
                % at each gap check pass save something O(n)
                if options.useCache
                    cacheSizePerPoint = zeros(param.n, 1);
                    activeSetSizePerPoint = zeros(param.n, 1);
                    for i_object = 1:param.n
                        cacheSizePerPoint(i_object) = cache{i_object}.cache_size;
                        activeSetSizePerPoint(i_object) = sum( cache{i_object}.alpha_vec(cache{i_object}.cache_entries) > 0 );
                    end
                    progress.cacheSizePerPoint = [progress.cacheSizePerPoint; {cacheSizePerPoint}];
                    progress.activeSetSizePerPoint = [progress.activeSetSizePerPoint; {activeSetSizePerPoint}];
                end
                progress.trueGapsPerPoint = [progress.trueGapsPerPoint; {gap_vec_true}];
                progress.heuristicGapsPerPoint = [progress.heuristicGapsPerPoint; {gap_vec_heuristic}];
            end
            % end of logging current status
        end 
        
        % do not forget to update the heuristic gaps (when the logging is already done)
        gap_vec_heuristic = gap_vec_true;
        
        if gap_true <= options.gap_threshold
            fprintf('Duality gap below threshold at pass: %.3f -- stopping!\n', num_passes)
            fprintf('Exact gap: %g, gap_threshold: %g\n', gap_true, options.gap_threshold)
            break;
        else
            fprintf('Exact duality gap (fw) = %g at pass = %d, time = %f\n', gap_true, i_pass, toc(tStart));
        end
    end
    % time-budget exceeded?
    t_elapsed = toc(tStart);
    if (t_elapsed/60 > options.time_budget)
        fprintf('Time budget exceeded.\n');
        break;
    end
end

if options.output_model_before_batch_step && exist('model_before_batch_step', 'var')
    model = model_before_batch_step;
end

end % solver_BCFW_hybrid

function options = defaultOptions(n)
options = struct;

% lambda
options.lambda = 1/n;

% stopping parameters
options.gap_threshold = 0.1;
options.num_passes = 200; % max number of passes through data
options.time_budget = inf;

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

% flag saying whether to demand computation of the true gap at convergence
options.true_gap_when_converged = true;

% If you want to use the output dual variables (cache) the following flag must always be false, 
% otherwise solver_BCFW_hybrid.m will be outputting dual variables inconsistent with the primal.
% However, we observed that after the batch FW step the primal objective significantly goes up (the gap is still going down). 
% If you actually want to use the output with batch FW steps (e.g. measure test error) we recommend to set this flag to true.
% Does not have any effect when options.do_batch_step == false.
options.output_model_before_batch_step = false;

% initialize the value of global gap used for the global cache criterion
% inf - always fail
% 0 - always pass
options.init_gap_value_for_cache = 1; % value 1 is correct only if the loss is bounded by one, otherwise heuristics until the first batch pass

% the logging level, how many things to store?
% - 0: no logging at all
% - 1: save models after each pass, and some data statistics, should not cause huge memory and speed overheads
% - 2: save everything reasonable, including O(n) numbers per each gap check pass
% - 3: crazy logging, save everything possible per point at each pass; huge computational overheads, use this only for specific plots
options.logging_level = 1;

end % defaultOptions


