function compile_BCFW

% add all the require paths
setup_BCFW;

% init
rootPath = fileparts( mfilename('fullpath') );
curPath = pwd;

% graphCutMex_BoykovKolmogorov for binary_segmentation
cd(fullfile(rootPath, 'applications', 'binary_segmentation', 'helpers', 'graphCutMex_BoykovKolmogorov'));
build_graphCutMex;
cd(curPath);

% distance_transform for pose_estimation
cd(fullfile(rootPath, 'applications', 'pose_estimation'));
build_distance_transform;
cd(curPath);

% crfChain for test_chunking
cd(fullfile(rootPath, 'applications', 'text_chunking'));
build_crfChain;
cd(curPath);

end

