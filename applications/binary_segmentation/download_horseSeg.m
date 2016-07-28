function download_horseSeg( data_path, dataset_version )
% This function downloads the OCR dataset: https://pub.ist.ac.at/~akolesnikov/HDSeg/
% Input: (optional) data_path - path where to put the downloaded data (default: <package root path>/data/horseSeg)
%        (optional) dataset_size - version of the dataset to download: 'small', 'medium', 'large', 'original' (default: 'small')

if ~exist('data_path', 'var') || isempty(data_path)
    % the default path
    root_path = fileparts( mfilename('fullpath') );
    while ~exist( fullfile( root_path, 'setup_BCFW.m' ), 'file' )
        root_path = fileparts( root_path );
    end
    
    data_path = fullfile(root_path, 'data', 'horseSeg' );
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
if strcmpi(dataset_version, 'original')
    fprintf('Downloading annotations...\n')
    download_wget('https://pub.ist.ac.at/~akolesnikov/HDSeg/HDSeg.tar', data_path);
    
    fprintf('Downloading superpixels...\n')
    download_wget('https://pub.ist.ac.at/~akolesnikov/HDSeg/data.tar', data_path );
    
    fprintf('Uncompressing files...\n');
    untar(fullfile(data_path, 'HDSeg.tar'), data_path);
    untar(fullfile(data_path, 'data.tar'), data_path);
    
    fprintf('To download the original images follow instructions in %s', fullfile(data_path,'HDSeg', 'README'))
elseif strcmpi(dataset_version, 'small')
    fprintf('Downloading processed %s dataset...\n', dataset_version)
    download_wget([url_release,'/','horseSeg_small.mat'], data_path);
elseif strcmpi(dataset_version, 'medium')
    fprintf('Downloading processed %s dataset...\n', dataset_version)
    download_wget([url_release,'/','horseSeg_medium.mat'], data_path);
elseif strcmpi(dataset_version, 'large')
    fprintf('Downloading processed %s dataset...\n', dataset_version)
    download_wget([url_release,'/','horseSeg_large.zip'], data_path);
    fprintf('Uncompressing files...\n');
    unzip(fullfile(data_path, 'horseSeg_large.zip'), data_path);
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