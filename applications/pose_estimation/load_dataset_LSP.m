function [ patterns_train, labels_train, patterns_test, labels_test ] = load_dataset_LSP( data_name, data_path )
% load_dataset_LSP prepares the LSP dataset in the BCFW format

% usage [ patterns_train, labels_train, patterns_test, labels_test ] = load_dataset_LSP( data_name, data_path )

% [1] Articulated Pose Estimation by a Graphical Model with Image Dependent
% Pairwise Relations, Chen and Yuille

% Inputs :
%       data_name : name of the required dataset ('LSP_small' or
%       'LSP_full')
%       data_path : path to the dataset
%
% Outputs : patterns_train, labels_train - features and labels of the
% training set
%           patterns_test, labels_test - idem for the test set

fname = fullfile(data_path, [data_name, '.mat']);

if ~(exist(fname, 'file'))
   error('You need to download the dataset %s (see README)\n', data_name);
else
   load(fname);
end




















