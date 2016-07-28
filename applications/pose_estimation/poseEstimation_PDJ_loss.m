function loss = poseEstimation_PDJ_loss(param, Y_truth, Y_predict)
%poseEstimation_PDJ_loss computes the number of correctly detected joints
% with 0.1 precision threshold
%
% loss = poseEstimation_sqdist_loss(param, Y_truth, Y_predict)
%
%  The loss is the average square distance between the joints.
%
% It is important that it is consistent with the loss function used in the loss-augmented decoding function!

% scale is defined as distance between left shoulder (index 3) and right
% hip (index 21)
scale               = norm(Y_truth.joints(3,:)-Y_truth.joints(21,:));
% for each joint compute the distance renormalized by 10% of the distance
% between left shoulder and right hip
distGTPredRescaled  = sum( (1./(1*scale)*(Y_predict.joints - Y_truth.joints)).^2, 2 );
% take max with one
lossTab             = min(distGTPredRescaled, ones(size(distGTPredRescaled)));
% take average over joints (1-loss) strongly correlated to the PDJ
loss                = mean(lossTab); 
