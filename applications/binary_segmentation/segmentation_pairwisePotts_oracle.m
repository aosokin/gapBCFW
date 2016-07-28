function labels_struct = segmentation_pairwisePotts_oracle(param, model, X, Y)
%segmentation_pairwisePotts_oracle does the loss-augmented decoding on a given example (X, Y) using model.w as parameter
%
% The model consist of the unary potentials for label 1 (the potential for label 0 is its negation) and Potts pairwise potentials
%
% If Y is not given, then standard prediction is done (i.e. MAP decoding without the loss term).
%
% Energy function:
% E(Y) = sum_i ( [y_i = 1]*(f_i^T w^U) - [y_i = 0] (f_i^T w^U) ) 
%      + sum_ij [y_i ~= y_j]*( f_{ij}^T w^P ) 
% Here w^U and w^P are unary and pairwise feature vectors, respectively. abs() is the elementwise absolute value function

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

% set the unary potentials
potentials_unary = zeros(num_nodes, 2);

% construct the connectivity system
pairwiseTerms = nan( num_edges, 4 );
edge_ends = double( X.edge_structure' + 1 ); % convert from 0-indexing to 1-indexing
pairwiseTerms(:, 1 : 2) = edge_ends;

% construct the pairwise weights
weights_unary = model.w(1 : num_unary_features); % weights for unary potentials, w^U
weights_pairwise = model.w( num_unary_features+1 : num_unary_features+num_pairwise_features); % weights for pairwise potentials, w^P

unary_impact = double( weights_unary(:)' * X.features );

% unary potentials
potentials_unary(:, 2) = unary_impact(:);
potentials_unary(:, 1) = -unary_impact(:);

% Potts pairwise potentials
pottsWeights = double(weights_pairwise(:)') * X.features_pairwise;

isNonSubmodular = pottsWeights < 0;
if sum(isNonSubmodular) > 0
    fprintf('CAUTION! %d non-submodular edges detected, truncating them\n', sum(isNonSubmodular));
    pottsWeights(isNonSubmodular) = 0;
end
pairwiseTerms(:, 3) = pottsWeights(:);
pairwiseTerms(:, 4) = pottsWeights(:);

if exist('Y', 'var') && ~isempty(Y)
    % adjust the unary potentials with the loss
    node_weights = Y.superpixel_size;
    node_weights = node_weights / sum( node_weights );
    
    switch param.lossType
        case 'hamming'
            unaryUpdate = getLossAugmentedUpdate_hamming( Y.labels, node_weights );
        case 'hammingBalanced'
            unaryUpdate = getLossAugmentedUpdate_hammingBalanced( Y.labels, node_weights );
        otherwise
            error('Unknown loss!')
    end
    potentials_unary = potentials_unary - unaryUpdate;
    
    labels_struct = Y;
else
    labels_struct = struct;
end

% do the decoding
unaryTerms = potentials_unary(:, [2, 1]);
[score_value, labels] = graphCutMex(unaryTerms, pairwiseTerms );

% output the structure as labels
labels_struct.labels = labels;

% % check the correctness
% phi = segmentation_pairwisePotts_featuremap(param, X, labels_struct);
% S = model.w' * phi;
% if exist('Y', 'var') && ~isempty(Y)
%     switch param.lossType
%         case 'hamming'
%             L = segmentation_hamming_loss(param, Y, labels_struct);
%         case 'hammingBalanced'
%             L = segmentation_hammingBalanced_loss(param, Y, labels_struct);
%         otherwise
%             error('Unknown loss!')
%     end
% else
%      L = 0;
% end
% fprintf('Current error is %f, non-submodular edges: %f%%\n', abs( -S-L - score_value ), 100 * sum(isNonSubmodular) / num_edges);
% 
end
