function loss = segmentation_hammingBalanced_loss(param, Y_truth, Y_predict)
%segmentation_hammingBalanced_loss computes the Hamming loss (with class imbalance taken into account) between the two predictions
%
% loss = segmentation_hammingBalanced_loss(param, Y_truth, Y_predict)
%
% The loss is normalized to [0,1], the sizes of superpixels zre taken into the account
%
% It is important that it is consistent with the loss function used in the loss-augmented decoding function!

if numel(Y_truth.labels) ~= numel(Y_predict.labels)
    error('The labelings passed into the loss function are not compatible');
end

is_foreground = Y_truth.labels(:) == 1;

size_foreground = sum( is_foreground .* Y_truth.superpixel_size(:) );
size_background = sum( ~is_foreground .* Y_truth.superpixel_size(:) );

loss_foreground = sum( ( is_foreground & (Y_predict.labels(:) == 0)) .* Y_truth.superpixel_size(:) ) / size_foreground;
loss_background = sum( (~is_foreground & (Y_predict.labels(:) == 1)) .* Y_truth.superpixel_size(:) ) / size_background;

loss = ( loss_foreground + loss_background ) / 2;

% check for the case when one of the labels is not present in the ground truth
if all(~is_foreground)
    loss = loss_background;
elseif all(is_foreground)
    loss = loss_foreground;    
end
