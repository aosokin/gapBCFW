function [w,wreg,w0,noneg] = model2vec(model)
% [w,wreg,w0,nonneg] = model2vec(model)
w     = zeros(model.len,1);
w0    = zeros(model.len,1);
wreg  = ones(model.len,1);
noneg = uint32([]);

for x = model.bias
  j = x.i:x.i+numel(x.w)-1;
  w(j) = x.w;
end

for x = model.apps
  j = x.i:x.i+numel(x.w)-1;
  w(j) = x.w;
  % Enforce
  w0(j) = .001;
  noneg = [noneg uint32(j)];
end

for x = model.pdefs
  j = x.i:x.i+numel(x.w)-1;
  w(j) = x.w;
  % Enforce
  w0(j) = .001;
  noneg = [noneg uint32(j)];
end

for x = model.gaus
  j = x.i:x.i+numel(x.w)-1;
  w(j) = x.w;
  % Enforce minimum quadratic deformation cost of 0.01
  j = [j(1),j(3)];
  w0(j) = .001;
  noneg = [noneg uint32(j)];
end
% Regularize root biases differently
for i = 1:length(model.components)
  b = model.components{i}(1).biasid;
  x = model.bias(b);
  j = x.i:x.i+numel(x.w)-1;
  wreg(j) = .001;
end
