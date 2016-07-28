function phi = poseEstimation_featuremap( param, X, Y )
% poseEstimation_featuremap computes the joint feature map phi(X, Y)
% 
% This code is inspired from the code of [1]:
%
% [1] Articulated Pose Estimation by a Graphical Model with Image Dependent
% Pairwise Relations, Chen and Yuille
%
% phi = poseEstimation_featuremap(param, X, Y)
% 
%
% param contains the information about the CNN model and also the graph
% structure of the tree
%
% X contains information about the sample :
% - X.iminfo.im = path to img
% - X.iminfo.isflip = is the image flipped (consider for data augmentation)
% - X.iminfo.r_degree = degree of rotation applied (consider for data augmentation)
% - X.iminfo.scale_x (scale information about image)
% - X.iminfo.scale_y (scale information about image)
% - X.path_feat = path to the features
% - X.name_feat = name of the features (in case user specify another path)

% Y is as struct containing the information on the label we are considering
% Y.joints  = coordinates of the joints in the original image
% Y.mix_id  = index of the joint relation (see t_ij in [1])

% get labels (joints and mixture ids)
labels  = Y;

% get parameters of the feature model
model_feat     = param.model_feat; % contains information about feature computation
components     = param.components;

% use only the first component
parts = components{1};

% Load some precomputed features
load(fullfile(param.feature_path, X.name_feat));

% initiate the joint feature vector
phi = zeros(model_feat.len, 1);

% feature component for the bias
phi(1) = 1; % same as them could put it bigger to not regularize it

% dimension of the models
n_po = numel(parts);

% scale used (only one in this work)
scale = pyra(1).scale;

% Start at the root
p      = 1;
x1     = (labels.joints(p,1)-1)/scale+1;
y1     = (labels.joints(p,2)-1)/scale+1;
phi(parts(p).appI) = unary_map{1}{p}(round(y1), round(x1));

[Ny, Nx] = size(unary_map{1}{p});

for p=2:n_po   
    % UNARY POTENTIAL 
    % get the coordinate inside the feature map 
    x1     = (labels.joints(p,1)-1)/scale+1;
    y1     = (labels.joints(p,2)-1)/scale+1;
    phi(parts(p).appI) = unary_map{1}{p}(min(Ny,max(1,round(y1))), min(Nx,max(1,round(x1))));
    
    
    % BINARY POTENTIAL
    % coordinate of child
    xc  = labels.joints(p,1);
    yc  = labels.joints(p,2);
    
    % parent of part p
    par = parts(p).parent;
    xp  = labels.joints(par,1);
    yp  = labels.joints(par,2);
    
    cbid = find(parts(p).nbh_IDs == parts(par).pid);
    pbid = find(parts(par).nbh_IDs == parts(p).pid);
    
    % a) appearance terms
    x2     = (labels.joints(par,1)-1)/scale+1;
    y2     = (labels.joints(par,2)-1)/scale+1;
    
    % child
    phi(parts(p).pdefI(cbid)) = ...
           idpr_map{1}{p}{cbid}(min(Ny,max(1,round(y1))), min(Nx,max(1,round(x1))),labels.mix_id{p}(cbid));
    
    % parent
    phi(parts(par).pdefI(pbid)) = ...
           idpr_map{1}{par}{pbid}(min(Ny,max(1,round(y2))), min(Nx,max(1,round(x2))),labels.mix_id{par}(pbid));
    
    % b) deformations terms
    xcn = round(x1);
    xpn = round(x2);
    ycn = round(y1);
    ypn = round(y2);
    
    % child-parent
    phi(parts(p).gauI{cbid}(labels.mix_id{p}(cbid)):parts(p).gauI{cbid}(labels.mix_id{p}(cbid))+3) = ...
        defvector(parts(p), xcn, ycn, xpn, ypn, labels.mix_id{p}(cbid), cbid);
    
    % parent-child
    phi(parts(par).gauI{pbid}(labels.mix_id{par}(pbid)):parts(par).gauI{pbid}(labels.mix_id{par}(pbid))+3) = ...
        defvector(parts(par), xpn, ypn, xcn, ycn, labels.mix_id{par}(pbid), pbid);
    
end


end

