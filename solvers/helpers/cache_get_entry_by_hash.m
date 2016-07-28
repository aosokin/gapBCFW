function [ cache, cache_index ] = cache_get_entry_by_hash( cache, yhash_i, psi_i, loss_i, ...
    cache_scores, max_cache_size, removeFromActiveSetIfFullCache )
% [ cache, cache_index ] = cache_get_entry_by_hash( cache, yhash_i, psi_i, loss_i, ...
%     cache_scores, max_cache_size, removeFromActiveSetIfFullCache )
%
% cache_get_entry_by_hash extracts an entry from the cache given the hash-key
% if the suitable entry is not found, the function updates the cache
%
% Inputs max_cache_size, removeFromActiveSetIfFullCache are optional

if ~exist('max_cache_size', 'var') || isempty(max_cache_size)
    max_cache_size = cache.cache_size;
end
if ~exist('removeFromActiveSetIfFullCache', 'var') || isempty(removeFromActiveSetIfFullCache)
    removeFromActiveSetIfFullCache = false;
end

% look for the correct coordinate in active set:
if cache.cached_labels.isKey(yhash_i)
    % the item is in the cache, no need to add it
    cache_index = cache.cached_labels(yhash_i);
else
    % check if there is space in cache or need to remove some corners
    if cache.cache_size < max_cache_size
        % can easily add element to cache
        active_index = cache.cache_size + 1;
        cache.cache_size = cache.cache_size + 1;
        cache.cache_entries = [cache.cache_entries; active_index];
    else
        % there is no space in cache, need either to kick someone out or to increase size of the cache
        % first, try to kick out non-active elementa
        removed_alpha = min( cache.alpha_vec(cache.cache_entries) );
        if removed_alpha == 0
            % cache contains elements that are not in the active set
            % pick for removal the non active element having the minimal score
            non_active_elements = find( cache.alpha_vec(cache.cache_entries) == 0 );
            if exist('cache_scores', 'var') && ~isempty(cache_scores)
                [~, bad_cache_element] = min( cache_scores( cache.cache_entries(non_active_elements) ) );
            else
                bad_cache_element = 1;
            end
            active_index = cache.cache_entries(non_active_elements(bad_cache_element));
        else
            % have to remove item from the active set.
            % have a choice: remove something from the active set or increase the cache size
            if ~removeFromActiveSetIfFullCache
                % decide to violate the cache size constraint
                active_index = cache.cache_size+1;
                cache.cache_size = cache.cache_size+1;
                cache.cache_entries = [cache.cache_entries; active_index];
                
                if cache.cache_size > numel( cache.alpha_vec )
                    % increase cache size by 20%: memory reallocation
                    extra_cache_size = max( round( 0.2 * cache.cache_size ), 1);
                    cache.l_vec = [ cache.l_vec; nan(extra_cache_size, 1) ];
                    if issparse( cache.AFeat )
                        cache.AFeat = [ cache.AFeat, sparse( size(cache.AFeat,1), extra_cache_size ) ];
                    elseif isa(cache.AFeat, 'single')
                        cache.AFeat = [ cache.AFeat, nan( size(cache.AFeat,1), extra_cache_size, 'single' ) ];
                    else
                        cache.AFeat = [ cache.AFeat, nan( size(cache.AFeat,1), extra_cache_size, 'double' ) ];
                    end
                    cache.alpha_vec = [ cache.alpha_vec; nan(extra_cache_size, 1) ];
                    cache.keys = [cache.keys; cell(extra_cache_size, 1)];
                end
            else
                % have to remove an item from an active set; pick the element with minimal cache score
                if exist('cache_scores', 'var') && ~isempty(cache_scores)
                    [~, bad_cache_element] = min( cache_scores(cache.cache_entries) );
                else
                    [~, bad_cache_element] = min( cache.alpha_vec(cache.cache_entries) );
                end
                active_index = cache.cache_entries( bad_cache_element );
                % raise a flag that the cache is not fully consistent: sum of alphas does not equal one
                cache.alpha_sum_to_one = false;
            end
        end
    end
    
    % remove the old elements from cache
    if ~isempty( cache.keys{active_index} )
        cache.cached_labels.remove( cache.keys{active_index} );
    end
    
    % adding FW corner to cache and the active set:
    cache.l_vec(active_index,1) = loss_i;
    cache.AFeat(:, active_index) = psi_i(:);
    cache.alpha_vec(active_index,1) = 0;
    cache.cached_labels(yhash_i) = active_index;
    cache.keys{active_index} = yhash_i;
    cache_index = cache.cache_entries(active_index);
end

end

