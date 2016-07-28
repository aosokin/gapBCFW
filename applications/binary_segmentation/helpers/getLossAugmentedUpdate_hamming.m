function unaryUpdate = getLossAugmentedUpdate_hamming( labels, node_weights )
%getLossAugmentedUpdate_hamming computes the update of the unary potentials 
% caused by the augmentation of the score function with the Hamming loss
%
% unaryUpdate = getLossAugmentedUpdate_hamming( labels, node_weights )
%
% Input:
%       labels - the ground-truth labeling, double[num_nodes, 1] of {0,1}
%       node_weights - the weights of the nodes based e.g. on the superpixel area, double[num_nodes, 1] > 0
%
% Output:
%       unaryUpdate - additive update of the unary potentials, double[num_nodes, 2]

K = 2;
N = numel( labels );

% get the update for the potentials
update_forLabel = 1 - eye(K);
unaryUpdate = update_forLabel(labels(:) + 1, :);
unaryUpdate = bsxfun(@times, unaryUpdate, node_weights(:));

if size( unaryUpdate, 1) ~= N && size( unaryUpdate, 2) ~= K
    error('Loss-augmented update for the Hamming loss produces inconsistent result');
end

end

