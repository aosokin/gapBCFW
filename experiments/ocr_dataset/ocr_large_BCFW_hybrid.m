function ocr_large_BCFW_hybrid(dataPath, resultPath, lambda, gap_threshold, num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, maxCacheSize, stepType, gap_check, rand_seed )


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
% We support two different settings for the dataset:
%   -'ocr' (OCR-small in ICML 2016 paper): only 1 of 10 folds in the training set
%   -'ocr2' (OCR-large in ICML 2016 paper): 9 of 10 folds in the training set
data_name = 'ocr2';  % OCR small dataset
[patterns_train, labels_train, patterns_test, labels_test] = loadOCRData(data_name, dataPath );

% create problem structure:
param = [];
param.patterns = patterns_train;
param.labels = labels_train;
param.lossFn = @chain_loss;
param.oracleFn = @chain_oracle;
param.featureFn = @chain_featuremap;
param.hashFn = @hash_ocr_sequence;
param.n = length( param.patterns );
phi1 = param.featureFn(param, param.patterns{1}, param.labels{1}); % use first example to determine dimension
param.d = length(phi1); % dimension of feature mapping
param.using_sparse_features = issparse(phi1);
param.cache_type_single =  true; % make cache smaller
clear phi1;


%% main parameters
setupName = ['BCFW_hybrid', ...
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

% choose in which mode to run the combined solver
% 0 - all steps
% 1 - only FW
% 2 - pairwise
% 3 - FW and away
options.stepType = stepType;

% other parameters
options.gap_check = gap_check; % how often to compute the true gap
options.rand_seed = rand_seed; % random seed

options.quit_passes_heuristic_gap = true; % quit the block-coordinate passes when the heuristic gap is below options.gap_threshold
options.true_gap_when_converged = true; % require true gap when method converges


% run the solver
if ~exist(modelFile, 'file')
    [model, ~, ~, progress] = solver_BCFW_hybrid(param, options);
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
