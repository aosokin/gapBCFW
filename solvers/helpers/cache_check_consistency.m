function is_consistent = cache_check_consistency( cache )
% is_consistent = cache_check_consistency( cache )
%
% cache_check_consistency checks cache data structure of one object for consistency.
% This functions performs partial check which is relatively fast. 
% If you are suspecting a cache bug try using cache_check_consistency_full.m instead.

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

% check cache.cached_labels
if (length( cache.cached_labels )~=cache.cache_size)
    fprintf('CAUTION (%s)! Consistency of cache is broken, number of entries is cached_labels does not equal cache_size\n', mfilename);
    is_consistent = false;
end

% check values of the dual variables
if any( cache.alpha_vec(cache.cache_entries) < -eps )
    fprintf('CAUTION (%s)! Consistency of cache is broken, some values of alpha_vec are negative\n', mfilename);
    is_consistent = false;    
end
if  any( cache.alpha_vec(cache.cache_entries) > 1+100*eps )
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

