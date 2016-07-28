function loss_augmented_scores = cache_get_scores( cache, w, cache_ids )
% loss_augmented_scores = cache_get_scores( cache, w, cache_ids )
%
% cache_get_scores computes the values H_i(y, w) over all the entries in the cache.
% The cache oracle consists in maximizing over these entries.

if ~exist('cache_ids', 'var') || isempty(cache_ids)
    cache_entries = cache.cache_entries(:);
else
    cache_entries = cache.cache_entries(cache_ids(:));
end

% cache.AFeat can be single, if this is the case convert it to double to do operations in double
loss_augmented_scores = cache.l_vec(cache_entries) - ( w' * double(cache.AFeat(:, cache_entries)) )';
loss_augmented_scores = loss_augmented_scores(:);

end

