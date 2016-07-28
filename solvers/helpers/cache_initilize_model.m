function [ cache, model ] = cache_initilize_model( param, max_cache_size, cache, model, lambda )
% [ cache, model ] = cache_initilize_model( param, max_cache_size, cache, model, lambda )
%
% cache_initilize_model creates model and cache data structures.
% The function either creates a new model with default initialization, reconstructs the model from cache,
% or re-initializes the model from both cache and the old model.

using_sparse_features = isfield(param, 'using_sparse_features') && param.using_sparse_features;
cache_type_single = isfield(param, 'cache_type_single') && param.cache_type_single;

if ~exist('cache', 'var') || isempty(cache) % standard initialization
    % cache: cell array, one for each example, containing the data-structure for the
    %        dual variables:
    %
    % cached_labels: this is a dictionary (containers.Map) mapping labeling to the index in the cache
    %       each entry of the dictionary is: (y_key; index_for_feature); y_key is computed using param.hashFn(ylabel)
    % AFeat: this is a size of cache x d matrix, each row storing the feature vector for a specific output
    % alpha_vec: each row is entry for alpha variable; 0 when inactive...
    % l_vec: each row contains the loss for this output
    % keys: each row gives the key of the corresponding output in the active set...
    % cache_size: size of the current cache
    % cache_entries: cache positions, occupied by some elements
    % alpha_sum_to_one: flag showing whether model thinks that alphas is cache sum to one, i.e. full active set is in the cache
    cache = cell(param.n,1);
    for i_object = 1:param.n
        ylabel = param.labels{i_object};
        % psi_i(y) := phi(x_i,y_i) - phi(x_i, y)

        cache{i_object} = struct;
        if using_sparse_features
            psi_i = sparse(param.d, 1);
            cache{i_object}.AFeat        = sparse( numel(psi_i), max_cache_size ); % preallocate memory for cache
            cache{i_object}.AFeat(:, 1)  = psi_i(:); % feature of ground truth (0!)
        elseif cache_type_single
            % economy mode: create cache of type double 
            psi_i = zeros(param.d, 1, 'single');
            cache{i_object}.AFeat = nan( numel(psi_i), max_cache_size, 'single' ); % preallocate memory for cache
            cache{i_object}.AFeat(:, 1) = psi_i(:); % feature of ground truth (0!)
        else
            % normal mode: create cache of type double
            psi_i = zeros(param.d, 1);
            cache{i_object}.AFeat = nan( numel(psi_i), max_cache_size ); % preallocate memory for cache
            cache{i_object}.AFeat(:, 1) = psi_i(:); % feature of ground truth (0!)
        end
        cache{i_object}.alpha_vec    = nan( max_cache_size, 1 );
        cache{i_object}.alpha_vec(1) = 1.0;
        cache{i_object}.l_vec = nan( max_cache_size, 1 );
        cache{i_object}.l_vec(1) = param.lossFn(param, ylabel, ylabel); % should be 0 but just in case... ;)
        yhash = param.hashFn(ylabel);
        cache{i_object}.keys = cell( max_cache_size, 1 );
        cache{i_object}.keys{1} = yhash;
        cache{i_object}.cached_labels = containers.Map({yhash},{1});
        cache{i_object}.cache_size = 1;
        cache{i_object}.cache_entries = 1;
        cache{i_object}.alpha_sum_to_one = true;
        
        
        if ~cache_check_consistency( cache{i_object} )
            error('Cache after initialization does not pass the consistency check. Something went wrong!');
        end
    end
    % initialize the model with zeros
    model = initialize_model_zeros(param);
    
else % WARM START
    % reuse previous model
    
    % check consistency of cache
    n = param.n;
    d = param.d;
    for i_object = 1 : n
        if ~cache_check_consistency( cache{i_object} )
            fprintf('CAUTION! cache of object %d used for warm start in %s is not consistent. You should definitely debug this!\n', i_object, mfilename );
        end
        if cache_type_single && ~using_sparse_features
            cache{i_object}.AFeat = single( cache{i_object}.AFeat );
        end
    end
    
    % get the model
    if exist('model', 'var') && ~isempty(model)
        fprintf('Using warm_start model...\n')
        
        % some checks of the consistency of the model cab be done here
        if ~isstruct( model )
            error('model provided for initialization is not of type structure');
        end
        if ~isfield(model, 'w')
            error('Field <w> is not present in the model provided for initialization');
        end
        if ~isfield(model, 'ell')
            error('Field <ell> is not present in the model provided for initialization');
        end
        if numel(model.w) ~= param.d
            error( 'model.w used for warm start does not have the right dimension')
        end
        if ~isfield(model, 'wMat')
            error('Field <wMat> is not present in the model provided for initialization');
        end
        if ~isfield(model, 'ellMat')
            error('Field <ellMat> is not present in the model provided for initialization');
        end
        if numel(size(model.wMat)) ~= 2 || size(model.wMat, 1) ~= d || size(model.wMat, 2) ~= n
            error( 'model.wMat used for warm start does not have the right dimension')
        end
        if numel(model.ellMat) ~= n
            error( 'model.ellMat used for warm start does not have the right dimension')
        end
        if isfield(param, 'positivity')
            if ~isfield(model, 'v')
                error('Field <v> is not present in the model provided for initialization');
            end
            if numel(model.v) ~= param.d
                error( 'model.v used for warm start does not have the right dimension')
            end
        end
        
    else
        if ~exist('lambda','var') || isempty(lambda)
            error('When initializing from cache either model or lambda have to be provided!');
        end
        model = struct;
        if using_sparse_features
            model.wMat = sparse(d, n);
        else
            model.wMat = zeros(d, n);
        end
        model.ellMat = zeros(n,1);
        
        % reconstructing w_i = A_i*alpha_i & ell_i:
        for i_object = 1:n
            I = cache{i_object}.cache_entries;
            model.wMat(:,i_object) = double(cache{i_object}.AFeat(:,I))*cache{i_object}.alpha_vec(I)/(lambda*n);
            model.ellMat(i_object) = double(cache{i_object}.l_vec(I))'*cache{i_object}.alpha_vec(I)/n;
        end
        
        if isfield(param, 'positivity')
            % keep a non trucanted version of w (A\alpha = v)
            model.v = sum(model.wMat,2);
            % check for sparsity (v is kept full)
            if issparse(model.v)
                model.v = full(model.v);
            end
            % making sure that the quantities are consistent
            v_truncated = model.v;
            v_truncated(param.positivity==1) = subplus(model.v(param.positivity==1));
            
            model.w = v_truncated;
            model.ell = sum(model.ellMat);
        else
            model.w = sum(model.wMat,2);
            % check for sparsity (w is kept full)
            if issparse(model.w)
                model.w = full(model.w);
            end
            model.ell = sum(model.ellMat);
        end
    end
end

end
