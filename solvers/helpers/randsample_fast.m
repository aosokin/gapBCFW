function i = randsample_fast(weights)
%randsample_fast is a speeded up version of MATLAB's randsample.m built-in function.
% The speed up is mostly achieved by cutting of unnecessary functionality, that slows things down.
%
% Usage:
% i = randsample_fast(weights);
% (Equaivalent to calling i = randsample(numel(weights),1,true,weights);)
%

n = numel(weights);
if n == 0
    error([mfilename, ': empty weight vector']);
end

sum_weights = sum(weights);
if ~(sum_weights > 0) || ~all(weights>=0) % catches NaNs
    error([mfilename, ': weights are not good']);
end
probs = weights(:)' / sum_weights;

        
edges = min([0 cumsum(probs)],1); % protect against accumulated round-off
edges(n+1) = 1; % get the upper edge exact
[~, i] = histc(rand, edges);

end