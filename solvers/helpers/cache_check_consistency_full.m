function is_consistent = cache_check_consistency_full( cache )
% is_consistent = cache_check_consistency_full( cache )
%
% cache_check_consistency_full checks cache data structure of one object for consistency.
% This functions performs full check but is quite slow. If you
% are not suspecting a cache bug use cache_check_consistency.m instead.

is_consistent = true;
% check consistency of allocated elements of cache       
num_allocated_elemets = numel( cache.alpha_vec );
if num_allocated_elemets ~= numel(cache.l_vec) || num_allocated_elemets ~= numel(cache.keys) || num_allocated_elemets ~= size(cache.AFeat,2)
    fprintf('CAUTION (%s)! Consistency of cache is broken, number of allocated elements is different for different components of the structure\n', mfilename);
    is_consistent = false;
end
if cache.cache_size > num_allocated_elemets
    fprintf('CAUTION (%s)! Consistency of cache is broken, cache_size larger than number of allocated elements\n', mfilename);
    is_consistent = false;
end

% check consistency of cache.cache_entries
if numel( cache.cache_entries ) ~= cache.cache_size
    fprintf('CAUTION (%s)! Consistency of cache is broken, number of cache_entries does not equal cache_size\n', mfilename);
    is_consistent = false;
end
if numel(unique(cache.cache_entries)) ~= cache.cache_size
    fprintf('CAUTION (%s)! Consistency of cache is broken, number of unique elements of cache_entries does not equal cache_size\n', mfilename);
    is_consistent = false;    
end
if ~isempty( setdiff( cache.cache_entries, 1 : num_allocated_elemets))
    fprintf('CAUTION (%s)! Consistency of cache is broken, elements of cache_entries are not in 1:num_allocated_elements\n', mfilename);
    is_consistent = false;    
end

% check that allocated elements have NaN and empty elements in the right places
non_cache_indices = setdiff( 1 : num_allocated_elemets, cache.cache_entries );
if any( isnan( cache.alpha_vec( cache.cache_entries ) ) ) || ~all(isnan(cache.alpha_vec(non_cache_indices)))
    fprintf('CAUTION (%s)! Consistency of cache is broken, alpha_vec contains NaN in the wrong places\n', mfilename);
    is_consistent = false;        
end
if any( isnan( cache.l_vec( cache.cache_entries ) ) ) || ~all(isnan(cache.l_vec(non_cache_indices)))
    fprintf('CAUTION (%s)! Consistency of cache is broken, l_vec contains NaN in the wrong places\n', mfilename);
    is_consistent = false;        
end

if ~issparse(cache.AFeat) && (any(any(isnan(cache.AFeat(:, cache.cache_entries)),1),2) || ~all(all(isnan(cache.AFeat(:, non_cache_indices)),1),2))
    fprintf('CAUTION (%s)! Consistency of cache is broken, AFeat contains NaN in the wrong places\n', mfilename);
    is_consistent = false;        
end
empty_keys = cellfun('isempty', cache.keys);
if any(empty_keys(cache.cache_entries)) || ~all(empty_keys(non_cache_indices))
    fprintf('CAUTION (%s)! Consistency of cache is broken, keys contain empty entries in the wrong places\n', mfilename);
    is_consistent = false;        
end

% check cache.cached_labels
if (length( cache.cached_labels )~=cache.cache_size)
    fprintf('CAUTION (%s)! Consistency of cache is broken, number of entries is cached_labels does not equal cache_size\n', mfilename);
    is_consistent = false;
end
all_values_dict = values(cache.cached_labels);
all_values_array = cat(1, all_values_dict{:});
if ~isempty(setxor(all_values_array, cache.cache_entries))
    fprintf('CAUTION (%s)! Consistency of cache is broken, values of cached_labels is not the same set as cache_entries\n', mfilename);
    is_consistent = false;
end

% check values of the dual variables
if any( cache.alpha_vec(cache.cache_entries) < -eps )
    fprintf('CAUTION (%s)! Consistency of cache is broken, some values of alpha_vec are negative\n', mfilename);
    is_consistent = false;    
end
if  any( cache.alpha_vec(cache.cache_entries) > 1+eps )
    fprintf('CAUTION (%s)! Consistency of cache is broken, some values of alpha_vec are larger than one\n', mfilename);
    is_consistent = false;    
end
if  any( sum(cache.alpha_vec(cache.cache_entries)) > 1+100*eps )
    fprintf('CAUTION (%s)! Consistency of cache is broken, sum of values of alpha_vec is larger than one\n', mfilename);
    is_consistent = false;    
end
if cache.alpha_sum_to_one && sum(cache.alpha_vec(cache.cache_entries)) < 1.0-1000*eps
    fprintf('CAUTION (%s): Consistency of cache is broken, alpha_sum_to_one flag is on and at the same time sum of alpha_vec is smaller that one by  %f\n', mfilename, 1 - sum(cache.alpha_vec(cache.cache_entries)) );
    is_consistent = false;    
end

end

