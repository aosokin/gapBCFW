function newResults = linearInterpolation( numPasses, curResults, numPassTarget )
%LINEARINTERPOLATION Summary of this function goes here
%   Detailed explanation goes here


newResults = curResults;
num_point_on_plots = numel( numPassTarget );
for i=1:numel(curResults)
    
    [val_sort,ind_sort] = sort(numPasses{i});
    results_for_plot = curResults{i}(ind_sort);
    
%     % do the median filtering if the number of point on the plots is too large
%     num_points_region_of_interest = sum( val_sort >= min(numPassTarget) & val_sort <= max(numPassTarget) );
%     medianRadius = ceil(num_points_region_of_interest/num_point_on_plots);
%     filteredSignal = medfilt1( results_for_plot, medianRadius );
%     
%     % do not do filtering at the borders to avoid edge effects
%     toChangeIndices = ceil(medianRadius/2) : 1 : numel(results_for_plot) - ceil(medianRadius/2);
%     if ~isempty( toChangeIndices )
%         toChangeIndices(toChangeIndices==0) = [];
%         results_for_plot(toChangeIndices) = filteredSignal(toChangeIndices);
%     end
    
    newResults{i} = interp1(val_sort,log(results_for_plot), numPassTarget);
    newResults{i} = exp(newResults{i});
end



end

