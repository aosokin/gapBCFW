function hLeg = plotMultipleMethodsWithErrorsBars( results, xTarget, nPtsMarker, alpha )
%PLOTMULTIPLEMETHODSWITHERRORSBARS Summary of this function goes here
%   Detailed explanation goes here
% results is a struct (one element per method)
% results.xValues   = cell array (for the different seeds)
% results.yValues   = cell array (for associated values)
% results.color     = color for the lines (RGB or 'r')
% results.marker    = marker style (will be filled by color)
% results.linestyle = style for the line
% results.name      = name for the legend
% xTarget is the target for x axis

% number of points that we draw
if ~exist('nPtsMarker', 'var')
    nPtsMarker = 4;
end

% alpha correspond to the quantile that we display for random seeds
if ~exist('alpha', 'var')
    alpha = 0;
end


nMethods     = numel(results);
methodNames  = cell(nMethods,1);
% plot the lines with error bars
hold on;
for m=1:nMethods
     methodNames{m} = results(m).name;
     
     % filter the results
     xValues = results(m).xValues;
     yValues = results(m).yValues;
     
     [xValue, yValue, errValue] = getFilterResults(xValues, yValues, xTarget, alpha); 
     % error bars
     if ~isequal(errValue(1,:),errValue(2,:))
          plotshaded(xValue, errValue, results(m).color);
     end
     % line plot
     plot(xValue, yValue, results(m).linestyle{:}, 'Color', results(m).color);    
end

% do the marker afterwards
plotIds = [];
for m=1:nMethods        
    % filter the results
    xValues = results(m).xValues;
    yValues = results(m).yValues;   
    [xValue, yValue, ~] = getFilterResults(xValues, yValues, xTarget, alpha); 
    
    % distribute them uniformly on numPass Target 
    % random starting point 
    stepUniform       = floor(numel(xValue)/(nPtsMarker+1));
    startSubSamplePts = floor(rand()*(stepUniform+1))-floor(stepUniform/2);
    subSamplePoints   = stepUniform:stepUniform:(numel(xValue)-stepUniform);
    subSamplePoints   = startSubSamplePts+subSamplePoints;

    % plot all markers without line
    plot(xValue(subSamplePoints), yValue(subSamplePoints), results(m).marker{:}, 'MarkerFaceColor', results(m).color, ...
        'Color',  results(m).color, 'LineStyle', 'none'); 

    % plot the last point to get the legend
    plotIds(m) = plot(0, 0,...
         results(m).linestyle{:},'MarkerFaceColor', results(m).color,results(m).marker{:},'Color', results(m).color); 
          
end

% display the legend
hLeg = legend(plotIds, methodNames);

end

