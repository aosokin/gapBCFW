% This demo applies the structured SVM to the LSP dataset collected in this paper:
%
% S. Johnson, M. Everingham
% Clustered Pose and Nonlinear Appearance Models for Human Pose Estimation 
% BMVC 2010
% Project webpage: http://www.comp.leeds.ac.uk/mat4saj/lsp.html
%
% The model used here is the one used in this paper:
%
% X. Chen, A. Yuille
% Articulated Pose Estimation by a Graphical Model with Image Dependent Pairwise Relations
% NIPS 2014
% Project webpage: 
% http://www.stat.ucla.edu/~xianjie.chen/projects/pose_estimation/pose_estimation.html
%
% The features needed for this script have been precomputed using the code
% of Chen and Yuille and can be downloaded here:
%
% TODO : put address of the features
%
% The model from Chen and Yuille is also required for this code, it can be
% downloaded here:
%
% TODO: put address of the CNN model
%
% The structured considered here is a tree corresponding to the human
% skeleton.
% The structured score is a function of binary variables containing unary 
% and pairwise potentials.
% The unary potentials are obtained from a CNN trained to classify
% individual joints.
% The binary potentials are composed of image dependent terms that come from a CNN
% and a geometric dependent terms.
% The max oracle for this problem is a max max-sum belief propagation 
% algorithm on an acyclic graph with messages computed with the generalized 
% distance transform (GDT)
%
% For the details of the model see Section I.4 of this paper:
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

data_path = fullfile(rootPath, 'data', 'LSP' );

%% load the dataset
% We support two different versions of the dataset: 'LSP_small', 'LSP_full'
%   - 'LSP_small' contains 100 images for training
%   - 'LSP_full' contains the original training set composed of 1000
%   images.
param           = struct;
param.data_path = data_path;
param.data_name = 'LSP_small'; % 'LSP_small' or 'LSP_full'

[param.patterns, param.labels, patterns_test, labels_test] = ...
    load_dataset_LSP(param.data_name, param.data_path);

%% setup the solver

% create problem structure:

% -----------------------
% LSP specific parameters
% -----------------------

% CNN model and feature path
model_feat           = load(fullfile(data_path, 'CNN_Deep_13_graphical_model.mat'));
param.model_feat     = model_feat.model;
[param.components,~] = modelcomponents(param.model_feat);
param.feature_path   = fullfile(data_path, 'features');

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

% -----------------------
% Training options
% -----------------------

options               = struct;
options.lambda        = 100;   % regularization parameter
options.gap_threshold = 0.01;  % duality gap stopping criterion
options.num_passes    = 100;     % max number of passes through data

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
