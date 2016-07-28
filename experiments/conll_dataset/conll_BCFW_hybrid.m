function conll_BCFW_hybrid(dataPath, resultPath, lambda, gap_threshold, num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, maxCacheSize, stepType, gap_check, rand_seed )

lambda
gap_threshold
num_passes
time_budget
sample
useCache
cacheNu
cacheFactor
maxCacheSize
stepType
gap_check
rand_seed

%% prepare the dataset
[patterns_train, labels_train, patterns_test, labels_test] = load_dataset_conll(dataPath);

%% create problem structure:
param              = [];
param.patterns     = patterns_train;
param.labels       = labels_train;
param.lossFn       = @chain_conll_loss;
param.oracleFn     = @chain_conll_oracle;
param.featureFn    = @chain_conll_featuremap;
param.hashFn       = @hash_conll_sequence;
param.n            = length(param.patterns);
phi1               = param.featureFn(param, param.patterns{1}, param.labels{1}); % use first example to determine dimension
param.d            = length(phi1); % dimension of feature mapping
param.using_sparse_features = issparse(phi1);


%% main parameters
setupName = ['BCFWHcombined', ...
    '_lambda',num2str(lambda), ...
    '_gapThreshold',num2str(gap_threshold), ...
    '_numPasses',num2str(num_passes),...
    '_timeBudget', num2str(time_budget), ...
    '_sample_',sample, ...
    '_useCache', num2str(useCache), ...
    '_cacheNu', num2str(cacheNu), ...
    '_cacheFactor', num2str(cacheFactor), ...
    '_maxCacheSize', num2str(maxCacheSize), ...
    '_stepType', num2str(stepType), ...
    '_gapCheck', num2str(gap_check), ...
    '_seed', num2str(rand_seed)];
    
%% prepare folder for the results
mkdir( resultPath );
resultFile = fullfile(resultPath, ['results_',setupName,'.mat']);
modelFile = fullfile(resultPath, ['models_',setupName,'.mat']);

%% run BCFW
options = struct;

% lambda
options.lambda = lambda;

% stopping parameters
options.gap_threshold = gap_threshold;
options.num_passes = num_passes; % max number of passes through data
options.time_budget = time_budget;

% gap sampling: 'uniform' or 'gap'
options.sample = sample;

% cache options
options.useCache = useCache;
options.cacheNu = cacheNu;
options.cacheFactor = cacheFactor;
options.maxCacheSize = maxCacheSize;

% kind of FW optimization 
options.stepType = stepType;

% other parameters
options.gap_check = gap_check; % how often to compute the true gap
options.rand_seed = rand_seed; % random seed

options.quit_passes_heuristic_gap = true; % quit the block-coordinate passes when the heuristic gap is below options.gap_threshold
options.true_gap_when_converged = true; % require true gap when method converges

% hack to save model only every gap_check iterations. For conll where the
% dimension is very high (1.6 million) but where we still keep the model as a non sparse
% vector (because w has no reason to be sparse a priori), we need to reduce
% the RAM used by doing this.
options.logging_level = 0.5; 

% run the solver
if ~exist(modelFile, 'file')
    [model, ~, ~, progress] = solver_BCFW_hybrid(param, options);
    
    % HACK for conll because we don't save all models
    
    % solve issues due to non saving all progress
    dim = size(progress.models{1}.w,1);
    % change things
    nPass = numel(progress.numPasses);
    % this will replace models
    models = cell(nPass, 1);
    verifiedGapModelId = progress.verifiedGapModelId;
    ct = 1;
    for j=gap_check:gap_check+1:nPass
        models{j} = progress.models{ct};
        verifiedGapModelId(ct) = j;
        ct = ct+1;
    end

    if ct==numel(verifiedGapModelId)
        verifiedGapModelId(end) = nPass-1;
        models{nPass-1} = progress.models{ct};
    end
    
    % add first model
    models{1}.w   = zeros(dim,1);
    models{1}.ell = 0;
    % add last model
    models{end} = progress.models{end};
    progress.models             = models;
    progress.verifiedGapModelId = verifiedGapModelId;
    save(modelFile, '-struct', 'progress', '-v7.3' );
else
    progress = load(modelFile);
end


% evaluate the models
if ~exist(resultFile, 'file')
    info = evaluate_models_fast( progress, param, options.lambda, param.patterns, param.labels, patterns_test, labels_test);
    info.lambda = progress.lambda;
    save(resultFile, '-struct', 'info', '-v7.3' );
end

