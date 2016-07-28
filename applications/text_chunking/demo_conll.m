% This demo applies the structured SVM to the CONLL dataset collected in this
% paper:
%
% Introduction to the CoNLL-2000 shared task: Chunking.
% Tjong Kim Sang, E. F. and Buchholz, S.
%
% The model used here is the one used in this paper:
%
% F. Sha and F.Pereira.
% Shallow parsing with conditional random fields
% NACL 2003
% Project webpage:
%
% More details on the graphical model used can be found at this webpage:
% http://www.chokkan.org/software/crfsuite/
%
% Otherwise see details of the model see Section I.2 of this paper:
%
% Anton Osokin, Jean-Baptiste Alayrac, Isabella Lukasewitz, Puneet K. Dokanian,
% Simon Lacoste-Julien
% Minding the Gaps for Block Frank-Wolfe Optimization of Structured SVMs
% ICML 2016
% Project page: http://www.di.ens.fr/sierra/research/gapBCFW

%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

data_path = fullfile(rootPath, 'data', 'CONLL' );

%% load the dataset
param           = struct;
[param.patterns, param.labels, patterns_test, labels_test] = ...
                                                  load_dataset_conll(data_path);

%% create problem structure:
param.lossFn       = @chain_conll_loss;
param.oracleFn     = @chain_conll_oracle;
param.featureFn    = @chain_conll_featuremap;
param.hashFn       = @hash_conll_sequence;
param.n            = length(param.patterns);
% use first example to determine dimension
phi1               = param.featureFn(param, param.patterns{1}, param.labels{1});
param.d            = length(phi1); % dimension of feature mapping
param.using_sparse_features = issparse(phi1);

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
    % standard prediction as don't give label as input
    ypredict = param.oracleFn(param, model, param.patterns{i});
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
