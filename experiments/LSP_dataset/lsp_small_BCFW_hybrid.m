function lsp_small_BCFW_hybrid(dataPath, resultPath, lambda, gap_threshold, num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, maxCacheSize, stepType, gap_check, rand_seed )


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
% We support two different versions of the dataset: 'LSP_small', 'LSP_full'
%   - 'LSP_small' contains 100 images for training (dataset used in
%   ICML2016)
%   - 'LSP_full' contains the original training set composed of 1000
%   images (dataset not used in ICML2016)

data_name = 'LSP_small';  % LSP small dataset
param = [];
[param.patterns, param.labels, patterns_test, labels_test] = ...
    load_dataset_LSP(data_name, dataPath);



param.data_path = dataPath;
param.data_name = data_name; 
% CNN model and feature path
model_feat           = load(fullfile(dataPath, 'CNN_Deep_13_graphical_model.mat'));
param.model_feat     = model_feat.model;
[param.components,~] = modelcomponents(param.model_feat);
param.feature_path   = fullfile(dataPath, 'features');
% positivity constraints
ind_gaus=[param.model_feat.gaus(:).i];
param.positivity = zeros(2677,1);

for j=1:numel(ind_gaus); 
    param.positivity(ind_gaus(j):2:ind_gaus(j)+3) = 1;
end

% -----------------
% Structured SVM 
% -----------------
% oracle, loss, feature (and hash function)
param.lossFn         = @poseEstimation_PDJ_loss;
param.oracleFn       = @poseEstimation_oracle;
param.featureFn      = @poseEstimation_featuremap;
param.hashFn         = @hash_pose_sequence;
% dimensions
param.n = length( param.patterns );
param.d = param.model_feat.len;
param.cache_type_single =  true; % make cache smaller


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
