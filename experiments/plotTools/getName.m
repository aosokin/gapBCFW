function name = getName(method, sampling, caching)

name = '%s, %s%s';

switch caching
    case 0
        cacheStr = '';
    case 1
        cacheStr = ', cache';
end

% choose in which mode to run the combined solver
% 0 - all steps
% 1 - only FW
% 2 - pairwise
% 3 - FW and away
switch method
    case 0
        name = sprintf(name, 'hybrid', sampling, cacheStr);
    case 1
        name = sprintf(name, 'BCFW', sampling, cacheStr);
    case 2
        name = sprintf(name, 'BCpFW', sampling, cacheStr);
    case 3
        name = sprintf(name, 'away', sampling, cacheStr);
end

