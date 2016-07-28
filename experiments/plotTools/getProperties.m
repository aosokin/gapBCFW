function [ color, marker, linestyle ] = getProperties( method, sampling, caching )
%GETPROPERTIES Summary of this function goes here
%   Detailed explanation goes here

switch caching
    case 0
        color = 'r';
        %color = [217,95,2]./255;
% <<<<<<< HEAD
%         color = [242,167,15]./255;
%     case 1
%         color = 'b';
%         %color = [67,162,202]./255;
%         color = [6,119,167]./255;
% =======
        %color = [242,167,15]./255;
        color = [224,163,42]./255;
    case 1
        color = 'b';
        %color = [67,162,202]./255;
        %color = [6,119,167]./255;
        color = [6,64,118]./255;
% >>>>>>> plotExp1
end

switch method
    case 1
        marker = {'Marker', 'square', 'MarkerEdgeColor', 'k'};
    case 2
        marker = {'Marker', 'o', 'MarkerEdgeColor', 'k'};
end

switch sampling
    case 'gap'
        linestyle = {'LineStyle', '-', 'LineWidth', 1};
    case 'uniform'
        linestyle = {'LineStyle', '--', 'LineWidth', 1.5};
end


end

