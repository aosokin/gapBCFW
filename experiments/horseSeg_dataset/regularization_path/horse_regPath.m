function horse_regPath(dataPath, resultPath, time_budget, ...
    true_gap, A, AplusB, stepType, eps_reg_path, gap_check, sample, useCache, ...
    num_passes, maxCacheSize, cacheNu, cacheFactor, rand_seed, datasetName)

dataPath
resultPath
time_budget
true_gap
A
AplusB
stepType
eps_reg_path
gap_check
sample
useCache
num_passes
maxCacheSize
cacheNu
cacheFactor
rand_seed
datasetName

%% prepare the dataset
param = struct;
param.data_path = dataPath;
param.data_name = datasetName; % 'horseSeg_small' or 'horseSeg_medium' or 'horseSeg_large'
if strcmpi(datasetName, 'horseSeg_large')
    [param.patterns, param.labels, patterns_test, labels_test] = load_dataset_horseSeg_featuresOnDisk(param.data_name, param.data_path);
else
    [param.patterns, param.labels, patterns_test, labels_test] = load_dataset_horseSeg(param.data_name, param.data_path);
end

param.oracleFn = @segmentation_pairwisePotts_oracle;
param.featureFn = @segmentation_pairwisePotts_featuremap;
param.hashFn = @hash_segmentation;
param.n = length( param.patterns );
phi1 = param.featureFn(param, param.patterns{1}, param.labels{1}); % use first example to determine dimension
param.d = length(phi1); % dimension of feature mapping
param.using_sparse_features = issparse(phi1);
param.cache_type_single =  true; % make cache smaller
clear phi1;

param.num_unary_features = 1969;
param.num_pairwise_features = 100;

param.positivity = zeros(param.d, 1);
param.positivity(param.num_unary_features + 1 : end) = 1;


param.lossType = 'hammingBalanced'; % 'hamming' or 'hammingBalanced'
switch param.lossType
    case 'hamming' 
        param.lossFn = @segmentation_hamming_loss;
    case 'hammingBalanced'
        param.lossFn = @segmentation_hammingBalanced_loss;
    otherwise
        error('Unknown loss!')
end

%% main parameters
setupName = ['regPath', ...
    '_numPasses',num2str(num_passes),...
    '_seed', num2str(rand_seed), ...
    '_gapCheck', num2str(gap_check), ...
    '_sample_',sample, ...
    '_useCache', num2str(useCache), ...
    '_cacheNu', num2str(cacheNu), ...
    '_cacheFactor', num2str(cacheFactor), ...
    '_stepType', num2str(stepType), ...
    '_timeBudget', num2str(time_budget), ...
    '_maxCacheSize', num2str(maxCacheSize), ...
    '_A', num2str(A), ...
    '_AplusB', num2str( AplusB ), ...
    '_epsRegPath', num2str(eps_reg_path),...
    '_trueGap', num2str(true_gap)];

%% prepare folder for the results
mkdir( resultPath );
resultFile = fullfile(resultPath, ['results_',setupName,'.mat']);
modelFile = fullfile(resultPath, ['models_',setupName,'.mat']);

%% run BCFW
% options structure:
options = struct;

% stopping parameters
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
options.true_gap_when_converged = true_gap; % require true gap when method converges

% multi-lambda parameters
options.quit_passes_heuristic_gap_eps_multiplyer = 0.8;

% no regularixation path
options.regularization_path = true;
options.regularization_path_min_lambda = 1e-5;
options.regularization_path_a = A;
options.regularization_path_b = AplusB - A;
if options.regularization_path_b <= 0
    fprintf( 'CAUTION (%s)! parameter be cannot be negative, A = %f, AplusB = %f!\n', mfilename, A, AplusB ) ;
    return;
end
options.regularization_path_eps = eps_reg_path;

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
