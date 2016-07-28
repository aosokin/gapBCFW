function cache = cache_set_dual_vars( cache, alphas )
% cache = cache_set_dual_vars( cache, alphas )
%
% cache_set_dual_vars is a setter function to update the alphas in the cache.
% MUST always call cache_get_dual_vars after cache_get_entry_by_hash prior to calling cache_set_dual_vars,
% bacause cache_get_entry_by_hash can update the cache.

if cache.alpha_sum_to_one && abs(sum(alphas) - 1) > 100*eps
    warning( ['CAUTION (%s)! after the update of dual variables they do not sum to one. Difference to one: ', mfilename, num2str(sum(alphas) - 1)] );
end

cache.alpha_vec( cache.cache_entries ) = alphas;

end

