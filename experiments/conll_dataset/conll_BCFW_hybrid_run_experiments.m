%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

dataPath = fullfile(rootPath, 'data', 'CONLL' );
resultPath = fullfile(rootPath, 'results' , 'conll_BCFW_hybrid' );

%% run the jobs
lambda        = {0.0001 ,0.01, 0.1};
useCache      = {1, 0};
gap_threshold = {10^-8};
num_passes    = {3000};
time_budget   = {10*60}; % 10 hours time limit
sample        = {'uniform', 'gap'};
cacheNu       = {0.01};
cacheFactor   = {0.25};
maxCacheSize  = {50};
stepType      = {1,2};
gap_check     = {10};
rand_seed     = num2cell(6:10);


task_launcher_combineArgs('conll_BCFW_hybrid', ...
    {dataPath},{resultPath}, lambda, gap_threshold, num_passes, time_budget,...
    sample, useCache, cacheNu, cacheFactor, ...
    maxCacheSize, stepType, gap_check, rand_seed);


% % Note that in our experiments, we used an APT toolbox designed for our cluster
% % (https://github.com/iXce/APT). Here is the exact code that we ran:
%
% APT_params();
% global APT_PARAMS;
% APT_PARAMS.exec_name = 'conll_BCFW_hybrid';
% APT_PARAMS.host_name = { 'mem001' };
% APT_compile( 'conll_BCFW_hybrid.m' );
% 
% APT_run('conll_BCFW_hybrid.m', ...
%     {dataPath},{resultPath}, lambdaValues, gapThresholdValues, numPassesValues, timeBudgetValues,...
%     samplingValues, cachingValues, cacheNu, cacheFactors, ...
%     maxCacheSizeValues, stepTypeValues, gapCheckValues, seedValues,...
%     'Memory', 6000, 'CombineArgs', 1, 'ClusterID', 2 );

