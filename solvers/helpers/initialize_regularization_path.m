function [lambda, model, gap_vec_heuristic, cache] = initialize_regularization_path(param, gap_threshold, construct_dual_vars, max_cache_size )
% [lambda, model, gap_vec_heuristic, cache] = initialize_regularization_path(param, gap_threshold, construct_dual_vars, max_cache_size )
%
% initialize_regularization_path implements the path initialization step described in Section 4 of ICML 2016 paper.

%% dataset details
d = param.d;
n = param.n;
phi = param.featureFn;
loss = param.lossFn;
oracle = param.oracleFn;

%% construct the zero solution: infinite lambda
if construct_dual_vars
    % create cache for zero-initialized model
    [cache, model] = cache_initilize_model( param, max_cache_size );
else
    cache = [];
    model = initialize_model_zeros(param);
end

%% find the maximum loss labelings
phi_max_loss = zeros(d,n);
phi_ground_truth = zeros(d,n);
y_i_max_loss = cell(n,1);
for i_object = 1:n
    % note, model should be initialized with zeros at this point, it respects positivity constrains
    y_i_max_loss{i_object} = oracle(param, model, param.patterns{i_object}, param.labels{i_object}); %max-oracle
    phi_max_loss(:, i_object) = phi(param, param.patterns{i_object}, y_i_max_loss{i_object});
    phi_ground_truth(:, i_object) = phi(param, param.patterns{i_object}, param.labels{i_object});
    model.ellMat(i_object) = loss(param, param.labels{i_object}, y_i_max_loss{i_object})/n;
end
% now we now model.w and model.wMat up to a scaling factor
psi_i = phi_ground_truth-phi_max_loss;
if isfield(param, 'cache_type_single') && param.cache_type_single
    % if cache is going to be single we need to truncate to single precision to avoid future in
    psi_i = single(psi_i);
    % need to cast back to double to have all the computations in double
    psi_i = double(psi_i);
end
model.wMat = psi_i/n;
model.ell = sum(model.ellMat);
if isfield(param, 'positivity')
    model.v = sum(model.wMat,2);
    % check for sparsity (v is kept full)
    if issparse(model.v)
        model.v = full(model.v);
    end
    model.w = model.v;
    model.w(param.positivity==1) = subplus(model.v(param.positivity==1));
else
    model.w = sum(model.wMat,2);
    % check for sparsity (w is kept full)
    if issparse(model.w)
        model.w = full(model.w);
    end
end

%% find labelings that correspond to parameters model.w
% max_y -w^T psiDiff(y) = max_y w^T phi(y) - w^T phi(y_i) = -w^T psiDiff(y_theta_i)
% where y_theta_i = argmax_y w^T phi(y)
theta = zeros(n,1); % theta(i_object) is going to be the maximum value of score with weight model.w
for i_object = 1:n
    y_i_theta = oracle(param, model, param.patterns{i_object});
    psi_i = phi_ground_truth(:,i_object) - phi(param, param.patterns{i_object}, y_i_theta);
    if isfield(param, 'cache_type_single') && param.cache_type_single
        % if cache is going to be single we need to truncate to single precision to avoid future in
        psi_i = single(psi_i);
        % need to cast back to double to have all the computations in double
        psi_i = double(psi_i);
    end
    theta(i_object) = -model.w'*psi_i;
end

%% compute lambda
lambda = (model.w'*model.w + sum(theta)/n) / gap_threshold;

%% assign primal variables
model.w = model.w / lambda;
model.wMat = model.wMat / lambda;
if isfield(param, 'positivity')
    model.v = model.v / lambda;
end

%% assign the gap upper bounds
gap_vec_heuristic = theta / (n*lambda) + (model.wMat)'*model.w *lambda;

%% assign dual variables
if  construct_dual_vars
    for i_object = 1:n
        yhash_i = param.hashFn( y_i_max_loss{i_object} );
        psi_i = phi_ground_truth(:, i_object) - phi_max_loss(:, i_object);
        ell_i = model.ellMat(i_object)*n;
        [ cache{i_object}, cache_index ] = cache_get_entry_by_hash( cache{i_object}, yhash_i, psi_i, ell_i, [], max_cache_size, false );
        alphas = cache_get_dual_vars( cache{i_object} );
        alphas = zeros(size(alphas));
        alphas(cache_index) = 1;
        cache{i_object} = cache_set_dual_vars( cache{i_object}, alphas );
    end
end

