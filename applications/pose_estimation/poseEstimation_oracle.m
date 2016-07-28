function labels_struct = poseEstimation_oracle(param, model, X, Y)
% segmentation_pairwisePotts_oracle does the loss-augmented decoding on a 
% given example (X, Y) using model.w as parameter
%
% [1] Articulated Pose Estimation by a Graphical Model with Image Dependent
% Pairwise Relations, Chen and Yuille
%
% This function is largely inspired from the detect function of [1].
% However it has been modified to our needs.
%
% If Y is not given, then standard prediction is done (i.e. MAP decoding
% without the loss term).

% labels_struct contains two terms :
% - labels_struct.joints = position of the joints in the original image
% - labels_struct.mix_id = id of the mixture components of the edges

% param contains the param of the CNN model (for feature computation) and
% also some info about the graphical model which is used (a tree)
% param.lossType = contains the type of the loss (here it's sqdist which is
% sum of square distance between joints which only modify unary potential)


% model.w contains the weights

% X contains information about the sample (taken from Chen and Yuile data
% structure)
% - X.iminfo.im       = path to img
% - X.iminfo.isflip   = is the image flipped (consider for data augmentation)
% - X.iminfo.r_degree = degree of rotation applied (consider for data augmentation)
% - X.iminfo.scale_x  = scale information about image (not used here)
% - X.iminfo.scale_y  = scale information about image (not used here)
% - X.name_feat       = name of the features (in case user specify another path)

% Y is as struct containing the information on the label we are considering
% Y.joints  = coordinates of the joints in the original image
% Y.mix_id  = index of the joint relation (see t_ij in [1])

% get parameters of the feature model
model_feat = param.model_feat; % contains information about feature computation
load(fullfile(param.feature_path, X.name_feat));

% scale considered. Note that in the original work [1], they considered
% several scale factor. Here we only keep one scale (the biggest one).
scale = pyra.scale;

% check if "loss augmented decoding" or not.
if exist('Y', 'var') && ~isempty(Y)
    % adjust the unary potentials with the loss
    unaryUpdate = getLossAugmentedUpdate_PDJ( Y, unary_map{1}, scale);
    labels_struct = Y;
else  
    labels_struct = struct;
end

% update the model with the weights (it updates the weights in the
% components to keep the same structure as [1]) 
model_feat         = vec2model(model.w, model_feat);
[components, apps] = modelcomponents(model_feat);

% start decoding
% use only the first component of the model
parts = components{1};
p_no  = numel(parts);

% Local scores
for p = 1:p_no
  % assign each deformation scores
  parts(p).defMap = idpr_map{1}{p};
  % --------------
  parts(p).appMap = unary_map{1}{p};
  f = parts(p).appid;
  if exist('Y', 'var') && ~isempty(Y)
        parts(p).score = parts(p).appMap * apps{f} + unaryUpdate{p};
  else
        parts(p).score = parts(p).appMap * apps{f};
  end
  parts(p).level = 1;
end

% Walk from leaves to root of tree, passing message to parent
for p = p_no:-1:2
  child  = parts(p);
  par    = parts(p).parent;
  parent = parts(par);
  cbid   = find(child.nbh_IDs == parent.pid);
  pbid   = find(parent.nbh_IDs == child.pid);
  
  % this part of the code uses distance transform (mex file)
  [msg,parts(p).Ix,parts(p).Iy,parts(p).Im{cbid},parts(par).Im{pbid}] ...
    = passmsg(child, parent, cbid, pbid);
  parts(par).score = parts(par).score + msg;
end

% Add bias to root score
parts(1).score = parts(1).score + parts(1).b;
rscore = parts(1).score;
[Y_m,X_m] = find(rscore >= max(rscore(:)));

%max(rscore(:))
% if ties choose the first one
breaktie = 1;
y = Y_m(breaktie);
x = X_m(breaktie);

[joints, mix_id] = backtrack_light(x,y,parts,scale);

labels_struct.joints = joints;
labels_struct.mix_id = mix_id;

% Backtrack through dynamic programming messages to estimate part locations
% and mixture id
function [joints, mix_id] = backtrack_light(x,y,parts,scale)
% slighty modified version of the backtrack from [1]
numparts = length(parts);
ptr      = zeros(numparts,2);
joints   = zeros(numparts,2);

mix_id = cell(numparts,1);

k   = 1;
ptr(k,:) = [x,y];
x1  = (x - 1)*scale+1;
y1  = (y - 1)*scale+1;

joints(k,:) = [x1 y1];

for k = 2:numparts
    p   = parts(k);
    par = p.parent;

    x   = ptr(par,1);
    y   = ptr(par,2);

    ptr(k,1) = p.Ix(y,x);
    ptr(k,2) = p.Iy(y,x);

    % update the joint position
    x1          = (ptr(k,1) - 1)*scale+1;
    y1          = (ptr(k,2) - 1)*scale+1;
    joints(k,:) = [x1 y1];

    cbid = find(p.nbh_IDs == parts(par).pid);
    pbid = find(parts(par).nbh_IDs == p.pid);

    cm = p.Im{cbid}(y,x);
    pm = parts(par).Im{pbid}(y,x);
    
    mix_id{k}(cbid)   = int32(cm);
    mix_id{par}(pbid) = int32(pm);
    
end


