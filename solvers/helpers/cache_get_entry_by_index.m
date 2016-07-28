function [ psi_cache, loss_cache, yhash_i ] = cache_get_entry_by_index( cache, cache_index )
% [ psi_cache, loss_cache, yhash_i ] = cache_get_entry_by_index( cache, cache_index )
%
% cache_get_entry_by_index is a getter function to extract elements from the cache given the index in the cache

cache_address = cache.cache_entries(cache_index);
psi_cache = double(cache.AFeat(:, cache_address)); % cache.AFeat can be single, if this is the case convert it to double to do operations in double
loss_cache = cache.l_vec(cache_address);
yhash_i = cache.keys{cache_address};

end

