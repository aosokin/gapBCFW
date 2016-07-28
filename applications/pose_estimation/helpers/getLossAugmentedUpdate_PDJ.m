function unaryUpdate = getLossAugmentedUpdate_PDJ( labels, unary_map, scale_cnn)
%getLossAugmentedUpdate_PDJ computes the update of the unary potentials 
% caused by the augmentation of the score function with the loss. The loss
% consists in a weighted sum of the squared distance between the ground truth
% joints and the predictions.
%
% unaryUpdate = getLossAugmentedUpdate_PDJ( labels, unary_map, scale_cnn )
%
% Input:
%       labels - the ground-truth labeling, struct with two main fields:
%                   
%                   * mix_id : ids of the pairwise joint
%                   * joints : 2D coordinates (in the image) of the joints
%                   (this is what is used to compute the loss)
%              - scale_cnn : images are down sampled when processed by the
%              CNN. This is the factor of down sampling which is used to
%              get the correspondences between feature map and image
%              coordinates.
%
% Output:
%       unaryUpdate - additive update of the unary potentials (n_nodes x 1 
%       cell)

n_po        = numel(unary_map);
unaryUpdate = unary_map;

% posMat contains the position (x,y) of each point in the feature map
[ly, lx] = size(unaryUpdate{1});
posMatX = repmat(1:lx, ly, 1);
posMatY = repmat((1:ly).', 1, lx);

% distance between the left shoulder and the right hip
scale_img  = norm(labels.joints(3,:)-labels.joints(21,:)); 

% transform the position in terms of coordinates in the original image rescaled by the
% normalized distance defined above

posMatX = ((posMatX-1)*scale_cnn+1)/(1*scale_img);
posMatY = ((posMatY-1)*scale_cnn+1)/(1*scale_img);

for p=1:numel(unary_map)
    xp = labels.joints(p,1)/(1*scale_img);
    yp = labels.joints(p,2)/(1*scale_img);
    unaryUpdate{p} = 1/(n_po)*min(((posMatX-xp).^2 + (posMatY-yp).^2), ones(size(posMatX)  )); 
end

end

