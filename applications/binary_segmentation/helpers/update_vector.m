function x = update_vector( x, index, update )
%update_vector adds numbers 'update' to positions 'index' of vector 'x'
% intuitivly this should be done as follows: 
% x(index) = x(index) + update; 
% but this operation does not work well when the index vector contains duplicates
%
% update can be either a vector of the same length as index or a scalar

updateSummed = accumarray(index, update, [length(x), 1], @sum);
x = x + updateSummed;

end
