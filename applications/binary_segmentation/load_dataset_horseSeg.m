function [patterns_train, labels_train, patterns_test, labels_test] = load_dataset_horseSeg(data_name, data_path)
%load_dataset_horseSeg prepares the horseSeg data set in the BCFW format
%
% [patterns_train, labels_train, patterns_test, labels_test] = load_dataset_horseSeg(data_name, data_path)
%
% Inputs:
%       data_name: name of the required dataset: 'horseSeg_small' or 'horseSeg_medium' or 'horseSeg_large' 
%       data_path: path to the downloaded data
%
% Outputs:
%       patterns_train, labels_train - features and labels of the training set    
%       patterns_test, labels_test - features and labels of the test set
%
% If you use this dataset please cite the following paper in any resulting publication:
%
% A. Kolesnikov, M. Guillaumin, V. Ferrari and C. H. Lampert 
% Closed-Form Approximate CRF Training for Scalable Image Segmentation 
% ECCV 2014
%
% The data can be download from here:
% https://pub.ist.ac.at/~akolesnikov/HDSeg/HDSeg.tar
% https://pub.ist.ac.at/~akolesnikov/HDSeg/data.tar

% load matlab-compatible dataset file, if it doesn't exist, create the .mat file
fname = fullfile(data_path, [data_name, '.mat']);
if (~exist(fname, 'file'))
    fprintf('Creating %s.mat for the first time...\n', data_name);
    [patterns_train, labels_train, patterns_test, labels_test] = convert_dataset_horseSeg( data_name, data_path );
    fprintf('Saving dataset to disk\n');
    save(fname, 'patterns_train', 'labels_train', 'patterns_test', 'labels_test', '-v7.3');
else
    fprintf('Loading already built %s dataset...\n', data_name)
    load(fname, 'patterns_train', 'labels_train', 'patterns_test', 'labels_test');
end
fprintf('Successfully loaded %s dataset.\n', data_name)

end

function [patterns_train, labels_train, patterns_test, labels_test] = convert_dataset_horseSeg( data_name, data_path )

%% setup the paths
config = struct;
config.trainFeatureFile = fullfile( 'data', 'horse-features.bin' );
config.testFeatureFile = fullfile( 'data', 'horse-features-test.bin' );
config.trainSuperPixelDir = fullfile( 'data', 'horse-superpixels' );
config.testSuperPixelDir = fullfile( 'data', 'horse-superpixels-test' );
config.imagePath = fullfile( 'HDSeg', 'images','HorseSeg' );
config.anotationPath = fullfile( 'HDSeg', 'annotations','HorseSeg' ); 

tStart = tic;

%% prepare the test set
fprintf('Converting the test data\n');
startIndex = 0;
endIndex = 241;
patterns_test = cell( endIndex - startIndex, 1 );
labels_test = cell( endIndex - startIndex, 1 );
               
curDataFile = fullfile( data_path, config.testFeatureFile );
% the following code saves the full description of hdf5 file into the tmp.txt
% description_hdf5 = evalc('h5disp( curDataFile )'); fid = fopen('tmp.txt', 'w'); fprintf(fid, '%s\n', description_hdf5); fclose(fid);

num_images = 0;
for iIndex = startIndex : 1 : endIndex - 1
    if mod(iIndex-startIndex+1, 100) == 0
        fprintf( 'Progress: %d of %d\n', iIndex-startIndex+1, endIndex - startIndex );
    end
    
    num_images = num_images + 1;
    indexStr = num2str( iIndex, '%06d');
    
    % get high-level information 
    labels_test{num_images} = struct;
    labels_test{num_images}.num_nodes = h5read(curDataFile, ['/superpixel_count_', indexStr] );
    labels_test{num_images}.num_states = 2;
    labels_test{num_images}.image_id = indexStr;
    
    patterns_test{num_images} = labels_test{num_images};
    
    % get the annotation
    curData = h5read(curDataFile, ['/annotation_file_', indexStr] );
    labels_test{num_images}.annotation_file = fullfile( config.anotationPath, 'test', curData{1} );
    if ~exist( fullfile( data_path, labels_test{num_images}.annotation_file ), 'file' )
        error( ['Cannot find the test annotation file: ', fullfile( data_path, labels_test{num_images}.annotation_file )] );
    end
    labels_test{num_images}.annotation_quality = 'manual';
    labels_test{num_images}.labels = h5read(curDataFile, ['/labels_', indexStr] );    
    
    % get the data
    patterns_test{num_images}.num_edges = h5read(curDataFile, ['/edge_count_', indexStr] );
    patterns_test{num_images}.edge_structure = h5read(curDataFile, ['/edge_structure_', indexStr] );
    patterns_test{num_images}.features = single( h5read(curDataFile, ['/features_', indexStr] ) );
    
    curData = h5read(curDataFile, ['/image_name_', indexStr] );
    [~, imageName, ~] = fileparts( curData{1} );
    imageFile = [imageName, '.JPEG'];
    imageFolder = 'test';
    patterns_test{num_images}.image_file = fullfile( config.imagePath, imageFolder, imageFile );
    if ~exist( fullfile( data_path, patterns_test{num_images}.image_file ), 'file' )
        error( ['Cannot find the test image file: ', fullfile( data_path, patterns_test{num_images}.image_file )] );
    end
    
    patterns_test{num_images}.superpixel_file = fullfile( config.testSuperPixelDir, [indexStr, '.png'] );
    if ~exist( fullfile( data_path, patterns_test{num_images}.superpixel_file), 'file' )
        error( ['Cannot find the test superpixel file: ', fullfile( data_path, patterns_test{num_images}.superpixel_file )] );
    end
    curSuperpixels = imread( fullfile( data_path, patterns_test{num_images}.superpixel_file ) );
    labels_test{num_images}.superpixel_size = accumarray( curSuperpixels(:) + 1, 1, [patterns_test{num_images}.num_nodes 1], @sum, 0);
    
    % check the data
    if size(patterns_test{num_images}.features, 2) ~= labels_test{num_images}.num_nodes
        error( ['Features of test image #', num2str(num_images), ' do not match its number of nodes'] );
    end
    if numel(labels_test{num_images}.labels) ~= labels_test{num_images}.num_nodes
        error( ['Labels of test image #', num2str(num_images), ' do not match its number of nodes'] );
    end
    if numel(labels_test{num_images}.superpixel_size) ~= labels_test{num_images}.num_nodes
        error( ['Superpixels of test image #', num2str(num_images), ' do not match its number of nodes'] );
    end
    
    % compute the pairwise features
    patterns_test{num_images}.features_pairwise = construct_pairwise_features( patterns_test{num_images}.features, patterns_test{num_images}.edge_structure+1 );
end

%% prepare the training set
fprintf('Converting the training data\n')
startIndex = 0;
switch data_name
    case 'horseSeg_small'
        endIndex = startIndex + 147;
    case 'horseSeg_medium'
        endIndex = startIndex + 147 + 5974;
    case 'horseSeg_large'
        endIndex = startIndex + 147 + 5974 + 19317;
    otherwise
        error(['Unknown dataset: ', data_name]);
end

patterns_train = cell( endIndex - startIndex, 1 );
labels_train = cell( endIndex - startIndex, 1 );

curDataFile = fullfile( data_path, config.trainFeatureFile );

num_images = 0;
for iIndex = startIndex : 1 : endIndex - 1
    if mod(iIndex-startIndex+1, 100) == 0
        fprintf( 'Progress: %d of %d\n', iIndex-startIndex+1, endIndex - startIndex );
    end
    
    num_images = num_images + 1;
    indexStr = num2str( iIndex, '%06d');
    
    % get high-level information 
    labels_train{num_images} = struct;
    labels_train{num_images}.num_nodes = h5read(curDataFile, ['/superpixel_count_', indexStr] );
    labels_train{num_images}.num_states = 2;
    labels_train{num_images}.image_id = indexStr;
    
    patterns_train{num_images} = labels_train{num_images};
    
    % get the annotation
    curData = h5read(curDataFile, ['/annotation_file_', indexStr] );
    [~, annotationName, ~] = fileparts( curData{1} );
    annotationParts = strsplit(annotationName, '_');
    % get the annotation type
    switch annotationParts{1}
        case 'annoa'
            labels_train{num_images}.annotation_quality = 'manual';
            annotationSuffix = '.gt';
        case 'annob'
            labels_train{num_images}.annotation_quality = 'bounding-box';
            annotationSuffix = '.bb';
        case 'annoc'
            labels_train{num_images}.annotation_quality = 'automatic';
            annotationSuffix = '.no';
        otherwise
            error( ['Unknown annotation type: ', annotationParts{1}]) ;
    end
    % get the annotation path
    annotationFile = [annotationParts{2}, '_', annotationParts{3}, annotationSuffix, '.png'];
    annotationFolder = annotationParts{2};
    if any( strcmpi(annotationFolder, {'2007', '2008', '2009', '2010', '2011'}) )
        annotationFolder = 'voc12';
    end
    labels_train{num_images}.annotation_file = fullfile( config.anotationPath, annotationFolder, annotationFile );
    if ~exist( fullfile( data_path, labels_train{num_images}.annotation_file), 'file' )
        error( ['Cannot find the train annotation file: ', fullfile( data_path, labels_train{num_images}.annotation_file )] );
    end
    
    labels_train{num_images}.labels = h5read(curDataFile, ['/labels_', indexStr] );    
    
    % get the data
    patterns_train{num_images}.num_edges = h5read(curDataFile, ['/edge_count_', indexStr] );
    patterns_train{num_images}.edge_structure = h5read(curDataFile, ['/edge_structure_', indexStr] );
    patterns_train{num_images}.features = single( h5read(curDataFile, ['/features_', indexStr] ) );
    
    curData = h5read(curDataFile, ['/image_name_', indexStr] );
    [~, imageName, ~] = fileparts( curData{1} );
    imageParts = strsplit(imageName, '_');
    if any( strcmp(imageParts{1}, {'imga', 'imgb', 'imgc'}) )
        if annotationParts{1}(end) ~= imageParts{1}(end)
            error( ['Image name does not match the annotation name: ', indexStr, '; ', imageName, '; ', annotationName] );
        end    
    else
        error( ['Unknown image type: ', imageParts{1}]) ;
    end
    if ~strcmp(imageParts{2}, annotationParts{2}) || ~strcmp(imageParts{3}, annotationParts{3})
        error(['Image and annotation do not match:',  imageParts{2}, '_', imageParts{3}, ' : ', annotationParts{2}, '_', annotationParts{3}]);
    end
    imageFile = [imageParts{2}, '_', imageParts{3}, '.JPEG'];
    imageFolder = imageParts{2};
    if any( strcmpi(imageFolder, {'2007', '2008', '2009', '2010', '2011'}) )
        imageFolder = 'voc12';
    end
    patterns_train{num_images}.image_file = fullfile( config.imagePath, imageFolder, imageFile );
    if ~exist( fullfile( data_path, patterns_train{num_images}.image_file), 'file' )
        error( ['Cannot find the train image file: ', fullfile( data_path, patterns_train{num_images}.image_file )] );
    end
    
    
    patterns_train{num_images}.superpixel_file = fullfile( config.trainSuperPixelDir, [indexStr, '.png'] );
    if ~exist( fullfile( data_path, patterns_train{num_images}.superpixel_file ), 'file' )
        error( ['Cannot find the train superpixel file: ', fullfile( data_path, patterns_train{num_images}.superpixel_file )] );
    end
    curSuperpixels = imread( fullfile( data_path, patterns_train{num_images}.superpixel_file ) );
    labels_train{num_images}.superpixel_size = accumarray( curSuperpixels(:) + 1, 1, [patterns_train{num_images}.num_nodes 1], @sum, 0);
    
    
	% check the data
    if size(patterns_train{num_images}.features, 2) ~= labels_train{num_images}.num_nodes
        error( ['Features of training image #', num2str(num_images), ' do not match its number of nodes'] );
    end
    if numel(labels_train{num_images}.labels) ~= labels_train{num_images}.num_nodes
        error( ['Labels of training image #', num2str(num_images), ' do not match its number of nodes'] );
    end
    if numel(labels_train{num_images}.superpixel_size) ~= labels_train{num_images}.num_nodes
        error( ['Superpixels of training image #', num2str(num_images), ' do not match its number of nodes'] );
    end
    
    % compute the pairwise features
    patterns_train{num_images}.features_pairwise = construct_pairwise_features( patterns_train{num_images}.features, patterns_train{num_images}.edge_structure+1 );
end
fprintf('Processing time: %fs\n', toc(tStart));

end

