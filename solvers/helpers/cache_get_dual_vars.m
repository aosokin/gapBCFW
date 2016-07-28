function alphas = cache_get_dual_vars( cache )
% alphas = cache_get_dual_vars( cache )
%
% cache_get_dual_vars is a getter function to extract alphas from cache

alphas = cache.alpha_vec( cache.cache_entries );

end

