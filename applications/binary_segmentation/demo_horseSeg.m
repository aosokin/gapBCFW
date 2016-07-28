% This demo applies the structured SVM to the HorseSeg dataset collected in this paper:
%
% A. Kolesnikov, M. Guillaumin, V. Ferrari, C. H. Lampert 
% Closed-Form Approximate CRF Training for Scalable Image Segmentation 
% ECCV 2014
% Project webpage: https://pub.ist.ac.at/~akolesnikov/HDSeg
%
% The data for this script can be downloaded from here:
%   https://pub.ist.ac.at/~akolesnikov/HDSeg/HDSeg.tar
%   https://pub.ist.ac.at/~akolesnikov/HDSeg/data.tar
%
% The structured score is a function of binary variables containing unary and pairwise potentials.
% The variables correspond to superpixels with label 0 corresponding to background and label 1
% corresponding to the foreground (horse) class.
% The unary potentials are formed from features representing superpixel and their neighborhoods.
% The pairwise potentials represent boundaries (pairwise comparison) between superpixels.
% The max oracle for this problem is a binary submodular energy minimization problem
% solved by the Boykov-Kolmogorov max-flow/min-cut algorithm.
% For the details of the model see Section I.3 of this paper:
%
% Anton Osokin, Jean-Baptiste Alayrac, Isabella Lukasewitz, Puneet K. Dokanian, Simon Lacoste-Julien
% Minding the Gaps for Block Frank-Wolfe Optimization of Structured SVMs
% ICML 2016
% Project page: http://www.di.ens.fr/sierra/research/gapBCFW

%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

data_path = fullfile(rootPath, 'data', 'horseSeg' );

%% load the dataset
% We support three different versions of the dataset: 'horseSeg_small', 'horseSeg_medium', and 'horseSeg_large'
%   - 'horseSeg_small' contains 147 images for training
%   - 'horseSeg_medium' contains 6,121 images for training
%   - 'horseSeg_large' contains 25,438 images for training
param = struct;
param.data_path = data_path;
param.data_name = 'horseSeg_small'; % 'horseSeg_small' or 'horseSeg_medium' or 'horseSeg_large'
[param.patterns, param.labels, patterns_test, labels_test] = load_dataset_horseSeg(param.data_name, param.data_path);
% % For large or even for medium version we recommend to use a version of the dataset that stores features
% % on disk and not in RAM. This significantly reduces the RAM requirements of the method, 
% % but reduces the oracle and feature computation time, because they need to read files from disk.
% [param.patterns, param.labels, patterns_test, labels_test] = load_dataset_horseSeg_featuresOnDisk(param.data_name, param.data_path);

%% setup the solver
% create problem structure
param.oracleFn = @segmentation_pairwisePotts_oracle;
param.featureFn = @segmentation_pairwisePotts_featuremap;
param.hashFn = @hash_segmentation;
param.n = length( param.patterns );
param.num_unary_features = 1969;
param.num_pairwise_features = 100;
param.d = param.num_unary_features + param.num_pairwise_features;
param.using_sparse_features = false;
param.cache_type_single =  true; % make cache smaller
param.positivity = zeros(param.d, 1);
param.positivity(param.num_unary_features + 1 : end) = 1; % add positivity contraint on the pairwise features

param.lossType = 'hammingBalanced'; % 'hamming' or 'hammingBalanced'
switch param.lossType
    case 'hamming' 
        param.lossFn = @segmentation_hamming_loss;
    case 'hammingBalanced'
        param.lossFn = @segmentation_hammingBalanced_loss;
    otherwise
        error('Unknown loss!')
end

% options structure
options = struct;
options.lambda = 100;
options.gap_threshold = 0.01; % duality gap stopping criterion
options.num_passes = 100; % max number of passes through data

%% run the solver
[model, ~, ~, progress] = solver_BCFW_hybrid( param, options );

%% loss on train set
avg_loss = 0;
for i=1:numel(param.patterns)
    ypredict = param.oracleFn(param, model, param.patterns{i}); % standard prediction as don't give label as input
    avg_loss = avg_loss + param.lossFn(param, param.labels{i}, ypredict);
end
avg_loss = avg_loss / numel(param.patterns);
fprintf('average loss on the training set: %f.\n', avg_loss);

% loss on test set
avg_loss = 0;
for i=1:numel(patterns_test)
    ypredict = param.oracleFn(param, model, patterns_test{i});
    avg_loss = avg_loss + param.lossFn(param, labels_test{i}, ypredict);
end
avg_loss = avg_loss / numel(patterns_test);
fprintf('average loss on the test set: %f.\n', avg_loss);
