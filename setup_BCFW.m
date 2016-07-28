function setup_BCFW

rootPath = fileparts( mfilename('fullpath') );

% add path to solvers
addpath( genpath( fullfile( rootPath, 'solvers' ) ) );

% add path the OCR code
addpath( genpath( fullfile( rootPath, 'applications', 'OCR' ) ) );

% add path to the binary segmentation code
addpath( genpath( fullfile( rootPath, 'applications', 'binary_segmentation' ) ) );

% add path to the text chunking
addpath( genpath( fullfile( rootPath, 'applications', 'text_chunking' ) ) );

% add path to the pose estimation
addpath( genpath( fullfile( rootPath, 'applications', 'pose_estimation' ) ) );

% add path to plotting tools
addpath( fullfile( rootPath, 'experiments', 'plotTools' ) );

end

