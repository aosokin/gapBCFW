function phi = segmentation_pairwisePotts_featuremap(param, X, Y)
%segmentation_pairwisePotts_featuremap computes the joint feature map phi(X, Y). [param is ignored]
%
% phi = segmentation_pairwisePotts_featuremap(param, X, Y)
% 
% It is important that the implementation is consistent with the segmentation_pairwisePotts_oracle.m!

% load feature from disk if not available
if ~isfield( X, 'features' ) || ~isfield(X, 'features_pairwise') || ~isfield(X, 'edge_structure')
    load( fullfile( param.data_path, X.feature_file ), 'features', 'features_pairwise', 'edge_structure' );
    X.features = features;
    X.features_pairwise = features_pairwise;
    X.edge_structure = edge_structure;
end

% problem dimensions
num_unary_features = size( X.features, 1 );
num_pairwise_features = size( X.features_pairwise, 1);
num_nodes = X.num_nodes;
num_states = X.num_states;
num_edges = X.num_edges;

if num_states ~= 2
    error('This oracle supports only binary problems!');
end

phi = zeros(1, num_unary_features + num_pairwise_features);

% get labels of the edges
labels = Y.labels(:);
edge_ends = double( X.edge_structure' + 1 ); % convert from 0-indexing to 1-indexing
label_first_end = labels( edge_ends(:, 1) );
label_second_end = labels( edge_ends(:, 2) );

% construct features for the unary potentials
phi(1 : num_unary_features) =  double( X.features * ((labels == 1) - (labels == 0)) );

% construct features for the pairwise potentials
phi( num_unary_features+1:num_unary_features+num_pairwise_features) = double(X.features_pairwise) * (label_first_end ~= label_second_end);

phi = -phi(:); % inverse the sign of features to make it consistent with the minimization oracle
if numel(phi) ~= num_unary_features + num_pairwise_features
    error( 'Incorrect size of the feature vector!' );
end
