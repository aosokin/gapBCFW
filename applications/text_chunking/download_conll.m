function download_conll( data_path )
% This function downloads the CONLL dataset preprocessed for BCFW structured
% SVM code. 
% Input: (optional) data_path - path where to put the downloaded data (default: <package root path>/data/LSP)

if ~exist('data_path', 'var') || isempty(data_path)
    % the default path
    root_path = fileparts( mfilename('fullpath') );
    while ~exist( fullfile( root_path, 'setup_BCFW.m' ), 'file' )
        root_path = fileparts( root_path );
    end
    
    data_path = fullfile(root_path, 'data', 'CONLL' );
end

% create the folder
if ~exist(data_path, 'dir')
    mkdir(data_path);
end

url_release = 'http://www.di.ens.fr/sierra/research/gapBCFW/release';

fprintf('Downloading processed CoNLL dataset...\n')
download_wget([url_release,'/','CONLL.zip'], data_path);
fprintf('Uncompressing files...\n');
unzip(fullfile(data_path, 'CONLL.zip'), data_path);

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