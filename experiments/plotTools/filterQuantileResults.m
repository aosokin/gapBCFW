function [ curveResults, shadedResults ] = filterQuantileResults( results, alpha )
%FILTERQUANTILERESULTS Summary of this function goes here
%   Detailed explanation goes here


results = results(:);
matValues = cell2mat(results);
nSample   = size(results,1);
nChop     = floor(nSample*alpha);

curveResults = median(matValues,1); % we plot the median
sortedValues = sort(matValues,1);
shadedResults = sortedValues([nChop+1, nSample-nChop],:);


end

