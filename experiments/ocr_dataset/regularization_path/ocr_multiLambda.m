function ocr_multiLambda(dataPath, resultPath, lambda_grid, gap_threshold, ...
    num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, ...
    maxCacheSize, stepType, gap_check, warm_start_type, rand_seed, datasetName )

lambda_grid
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
warm_start_type
rand_seed
datasetName


%% prepare the dataset
% OCR dataset:
% We support two different settings for the dataset (ocr: only one fold in
% training set, ocr2: all but one fold in training set
% -- ocr2 is the one that we have used in our experiments in the
% ICML 2013 paper)
if strcmpi(datasetName, 'large')
    data_name = 'ocr2';  % OCR large dataset
elseif strcmpi(datasetName, 'small')
    data_name = 'ocr';  % OCR small dataset
else
    error('Unknown dataset name.')
end
%[patterns_train, labels_train, patterns_test, labels_test] = loadOCRData(data_name, dataPath );
[patterns_train, labels_train] = loadOCRData(data_name, dataPath );

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
setupName = ['gridSearch', ...
    '_lambdaGrid_',lambda_grid, ...
    '_gapThreshold',num2str(gap_threshold), ...
    '_warmStart', warm_start_type, ...
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
lambda_params = sscanf(lambda_grid, 'factor%fexp%fto%f');
options.lambda_values = lambda_params(1).^(lambda_params(2):-1:lambda_params(3));

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

% multi-lambda parameters
options.quit_passes_heuristic_gap_eps_multiplyer = 0.8;
options.warm_start_type = warm_start_type; % 'keep_primal' or 'keep_dual' or 'none';

% no regularixation path
options.regularization_path = false;

% for debugging
options.check_lambda_change = false;

% run the solver
if ~exist(modelFile, 'file')
    [model, ~, ~, progress] = solver_multiLambda_BCFW_hybrid(param, options);
    save(modelFile, 'progress', '-v7.3' );
else
    load(modelFile, 'progress');
end

% evaluate the models
if ~exist(resultFile, 'file')
    info = evaluate_regPath_models( param, progress );
    save(resultFile, 'info', '-v7.3' );
end
