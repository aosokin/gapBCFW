function loss = segmentation_hamming_loss(param, Y_truth, Y_predict)
%segmentation_hamming_loss computes the Hamming loss between the two predictions
%
% loss = segmentation_hamming_loss(param, Y_truth, Y_predict)
%
% The loss is normalized to [0,1], the sizes of superpixels zre taken into the account
%
% It is important that it is consistent with the loss function used in the loss-augmented decoding function!

if numel(Y_truth.labels) ~= numel(Y_predict.labels)
    error('The labelings passed into the loss function are not compatible');
end

loss = sum( (Y_truth.labels(:) ~= Y_predict.labels(:)) .* Y_truth.superpixel_size(:) ) / sum( Y_truth.superpixel_size(:) );
