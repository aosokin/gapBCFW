% This demo applies the structured SVM to the OCR dataset by Ben Taskar. The structured
% model considered here is the standard chain graph, with the pixel values of
% the digit as unary features and a transition matrix of size num_states^2 as
% a pairwise potential. Additionally, we include a unary bias term for the first
% and last symbol in the sequence.

%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

dataPath = fullfile(rootPath, 'data', 'OCR' );

%% load the dataset
% We support two different settings for the dataset:
%   -'ocr' (OCR-small in ICML 2016 paper): only 1 of 10 folds in the training set
%   -'ocr2' (OCR-large in ICML 2016 paper): 9 of 10 folds in the training set
data_name = 'ocr'; %'ocr' or 'ocr2'
[patterns_train, labels_train, patterns_test, labels_test] = loadOCRData(data_name, dataPath);

%% setup the solver
% create problem structure
param = [];
param.patterns = patterns_train;
param.labels = labels_train;
param.lossFn = @chain_loss; %Hamming loss normalized to [0,1]
param.oracleFn = @chain_oracle; %max-product type decoding
param.featureFn = @chain_featuremap;
param.hashFn = @hash_ocr_sequence;
param.n = length( param.patterns );
phi1 = param.featureFn(param, param.patterns{1}, param.labels{1}); % use first example to determine dimension
param.d = length(phi1); % dimension of feature mapping
param.using_sparse_features = issparse(phi1);
clear phi1;

% options structure
options = struct;
options.lambda = 0.0373;
options.gap_threshold = 0.1; % duality gap stopping criterion
options.num_passes = 100; % max number of passes through data

%% run the solver
[model, ~, ~, progress] = solver_BCFW_hybrid( param, options );

%% loss on train set
avg_loss = 0;
for i=1:numel(patterns_train)
    ypredict = param.oracleFn(param, model, patterns_train{i}); % standard prediction as don't give label as input
    avg_loss = avg_loss + param.lossFn(param, labels_train{i}, ypredict);
end
avg_loss = avg_loss / numel(patterns_train);
fprintf('average loss on the training set: %f.\n', avg_loss);

% loss on test set
avg_loss = 0;
for i=1:numel(patterns_test)
    ypredict = param.oracleFn(param, model, patterns_test{i});
    avg_loss = avg_loss + param.lossFn(param, labels_test{i}, ypredict);
end
avg_loss = avg_loss / numel(patterns_test);
fprintf('average loss on the test set: %f.\n', avg_loss);

%% plot the progress of the solver
plot(progress.time(progress.verifiedGapModelId), progress.verifiedGap, 'b');
xlabel('time (s)');
ylabel('Frank-Wolfe duality gap');

