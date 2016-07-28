function [ convertTab ] = transform_time( tab, format )
%TRANSFORM_TIME Summary of this function goes here
%   Detailed explanation goes here

switch format
    case 'h'
        convertTab = tab./3600;
    case 'm'
        convertTab = tab./60;
    case 's'
        convertTab = tab;
end

end

