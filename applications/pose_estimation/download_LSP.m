function download_LSP( data_path, dataset_version )
% This function downloads the LSP dataset (http://www.comp.leeds.ac.uk/mat4saj/lsp.html)
% preprocessed for BCFW structured SVM code. The preprocessed features have
% been obtained with the code of Chen and Yuille
% (http://www.stat.ucla.edu/~xianjie.chen/projects/pose_estimation/pose_estimation.html)
% It will also download a part of the CNN model used by the code of Chen
% and Yuille.
% Input: (optional) data_path - path where to put the downloaded data (default: <package root path>/data/LSP)
%        (optional) dataset_size - version of the dataset to download: 'small', 'full' (default: 'small')

if ~exist('data_path', 'var') || isempty(data_path)
    % the default path
    root_path = fileparts( mfilename('fullpath') );
    while ~exist( fullfile( root_path, 'setup_BCFW.m' ), 'file' )
        root_path = fileparts( root_path );
    end
    
    data_path = fullfile(root_path, 'data', 'LSP' );
end
if ~exist('dataset_version', 'var') || isempty(dataset_version)
    dataset_version = 'small';
end


% create the folder
if ~exist(data_path, 'dir')
    mkdir(data_path);
end

url_release = 'http://www.di.ens.fr/sierra/research/gapBCFW/release';

% downloading
if strcmpi(dataset_version, 'small')
    fprintf('Downloading processed %s dataset...\n', dataset_version)
    download_wget([url_release,'/','LSP_small.mat'], data_path);
    
    if ~exist(fullfile(data_path,'CNN_Deep_13_graphical_model.mat'),'file')
        fprintf('Downloading model from Chen&Yuille...\n')
        download_wget([url_release,'/','CNN_Deep_13_graphical_model.mat'], data_path);
    end
    
    fprintf('Downloading individual features for %s...\n', dataset_version)
    download_wget([url_release,'/','features_LSP_small.zip'], data_path);
    fprintf('Uncompressing files...\n');
    unzip(fullfile(data_path, 'features_LSP_small.zip'), data_path);
    
elseif strcmpi(dataset_version, 'full')
    fprintf('Downloading processed %s dataset...\n', dataset_version)
    download_wget([url_release,'/','LSP_full.mat'], data_path);
    
    if ~exist(fullfile(data_path,'CNN_Deep_13_graphical_model.mat'),'file')
        fprintf('Downloading model from Chen&Yuille...\n')
        download_wget([url_release,'/','CNN_Deep_13_graphical_model.mat'], data_path);
    end
    
    fprintf('Downloading individual features for %s (big file ahead)...\n', dataset_version)
    download_wget([url_release,'/','features_LSP_full.zip'], data_path);
    fprintf('Uncompressing files...\n');
    unzip(fullfile(data_path, 'features_LSP_full.zip'), data_path);
else
    fprintf('Unknown dataset version: %s\n', dataset_version);
end

end


function exit_code = download_wget( URL, data_path )

wget_cmd = ['wget -P ', data_path, ' ', URL];
exit_code = system(wget_cmd);
if exit_code ~= 0
    fprintf('Downloading with MATLAB, no progress report...\n');
    [~, file_name, file_extension] = fileparts(URL);
    [~, status] = urlwrite(URL, fullfile(data_path, [file_name, file_extension]));
    exit_code = 1-status;
end

end