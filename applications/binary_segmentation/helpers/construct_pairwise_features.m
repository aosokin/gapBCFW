function features_pairwise = construct_pairwise_features( features_unary, edge_structure )

histogram_feature_indices = cell(9, 1);

% position features
histogram_feature_indices{1} = 1 + (1 : 16);
histogram_feature_indices{2} = 1 + 16 + (1 : 16);
histogram_feature_indices{3} = 1 + 32 + (1 : 16);

% color features
histogram_feature_indices{4} = 1 + 48 + (1 : 128);
histogram_feature_indices{5} = 1 + 48 + 128 + (1 : 128);
histogram_feature_indices{6} = 1 + 48 + 256 + (1 : 128);

% sift features
histogram_feature_indices{7} = 1 + 48 + 384 + (1 : 512);
histogram_feature_indices{8} = 1 + 48 + 384 + 512 + (1 : 512);
histogram_feature_indices{9} = 1 + 48 + 384 + 1024 + (1 : 512);

num_nodes = size( features_unary, 2 );
num_edges = size(edge_structure, 2);
num_hists = numel(histogram_feature_indices);

weight_possibilities = 2.^(-5 : 5);
num_weights = numel(weight_possibilities);
features_pairwise = nan( num_hists*num_weights, num_edges );
for i_hist = 1 : numel(histogram_feature_indices)
    X = features_unary(histogram_feature_indices{i_hist}, :).^2'; % chi2 distance have to be applied to L1 normalized hists, but ECCV 2014 paper has L2-normalized features
    distance_matrix = pdist2( X, X, 'chisq' );
    hist_distance = distance_matrix( edge_structure(1, :) + num_nodes * (edge_structure(2, :) - 1) );
   
    features_pairwise( (1 : num_weights) + (i_hist-1)*num_weights, :) = exp( -(hist_distance(:) * weight_possibilities)' );
end

% add a constant feature
features_pairwise = [ones(1, num_edges);  features_pairwise];

% convert feature to single to save space
features_pairwise = single(features_pairwise);

end

