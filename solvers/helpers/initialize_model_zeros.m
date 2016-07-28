function model = initialize_model_zeros(param)
% model = initialize_model_zeros(param)
%
% initialize_model_zeros creates a model data structure and initializes it with zeros (default initialization)

n = param.n; % number of training examples
d = param.d; % dimension of feature mapping

if isfield(param, 'using_sparse_features') && param.using_sparse_features
    model.wMat = sparse(d,n);
else
    model.wMat = zeros(d,n);
end

model.w = zeros(d,1);
model.ell = 0;
model.ellMat = zeros(n,1);

if isfield(param, 'positivity')
    % keep a non truncated version of w
    model.v = model.w;
end

end