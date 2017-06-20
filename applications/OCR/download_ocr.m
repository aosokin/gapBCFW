function download_ocr( data_path )
% This function downloads the OCR dataset: http://ai.stanford.edu/~btaskar/ocr/
% Input: (optional) data_path - path where to put the downloaded data (default: <package root path>/data/OCR)

if ~exist('data_path', 'var') || isempty(data_path)
    % the default path
    root_path = fileparts( mfilename('fullpath') );
    while ~exist( fullfile( root_path, 'setup_BCFW.m' ), 'file' )
        root_path = fileparts( root_path );
    end
    
    data_path = fullfile(root_path, 'data', 'OCR' );
end

% create the folder
if ~exist(data_path, 'dir')
    mkdir(data_path);
end

% downloading
download_wget('http://ai.stanford.edu/~btaskar/ocr/letter.names', data_path );
download_wget('http://ai.stanford.edu/~btaskar/ocr/letter.data.gz', data_path);

% uncompress
gunzip(fullfile(data_path, 'letter.data.gz'), data_path);

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
