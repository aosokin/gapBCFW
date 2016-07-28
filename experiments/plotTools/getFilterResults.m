function [ xValue, yValue, errValue ] = getFilterResults( xValues, yValues, xTarget, alpha )
%GETFILTERRESULTS Summary of this function goes here
%   Detailed explanation goes here

minVal = -inf;
maxVal = +inf;

for i=1:numel(xValues)
    minVal = max([min(xValues{i}), minVal]);
    maxVal = min([max(xValues{i}), maxVal]);
end

indStart = find(xTarget>=minVal,1,'first');
indEnd   = find(xTarget<=maxVal,1,'last');
xValue   = xTarget(indStart:indEnd);

newResults       = linearInterpolation(xValues, yValues, xValue);
[yValue, errValue] = filterQuantileResults(newResults, alpha);

 
end

