function unaryUpdate = getLossAugmentedUpdate_hammingBalanced( labels, node_weights )
%getLossAugmentedUpdate_hammingBalanced computes the update of the unary potentials 
% caused by the augmentation of the score function with the Hamming loss balanced over classes
%
% unaryUpdate = getLossAugmentedUpdate_hammingBalanced( labels, node_weights )
%
% Input:
%       labels - the ground-truth labeling, double[num_nodes, 1] of {0,1}
%       node_weights - the weights of the nodes based e.g. on the superpixel area, double[num_nodes, 1] > 0
%
% Output:
%       unaryUpdate - additive update of the unary potentials, double[num_nodes, 2]

K = 2;
N = numel( labels );

% get the weights for the nodes
is_foreground = labels == 1;
loss_max = sum( node_weights );
loss_foreground = sum( node_weights( is_foreground ) );
loss_background = sum( node_weights( ~is_foreground ) );

foreground_factor = loss_max / loss_foreground / 2;
background_factor = loss_max / loss_background / 2;
% check for the case if one of the labels is not present
if all( ~is_foreground )
    background_factor = background_factor * 2;
elseif all( is_foreground )
    foreground_factor = foreground_factor * 2;
end
    
node_weights_final = node_weights;
node_weights_final( is_foreground ) = node_weights_final( is_foreground )  * foreground_factor;
node_weights_final( ~is_foreground ) = node_weights_final( ~is_foreground ) * background_factor;

% get the update for the potentials
update_forLabel = 1 - eye(K);
unaryUpdate = update_forLabel(labels(:) + 1, :);
unaryUpdate = bsxfun(@times, unaryUpdate, node_weights_final);

if size( unaryUpdate, 1) ~= N && size( unaryUpdate, 2) ~= K
    error('Loss-augmented update for the Hamming loss produces inconsistent result');
end

end

