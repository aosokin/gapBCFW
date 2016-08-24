%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

dataPath = fullfile(rootPath, 'data', 'horseSeg' );
resultPath = fullfile(rootPath, 'results' , 'horseSeg_medium_BCFW_hybrid' );

%% run the job
lambda          = {10, 1, 1e-1};
gap_threshold   = {1e-6};
num_passes      = {20000};
time_budget     = {6 * 60};
sample          = {'gap', 'uniform'};
useCache        = {1, 0};
cacheNu         = {0.01};
cacheFactor     = {0.25};
maxCacheSize    = {50};
stepType        = {1, 2};
gap_check       = {10};
rand_seed       = num2cell(6 : 10);
datasetName     = {'horseSeg_medium'};

task_launcher_combineArgs('horse_BCFW_hybrid', ...
    {dataPath}, {resultPath}, lambda, gap_threshold, num_passes, time_budget, ...
    sample, useCache, cacheNu, cacheFactor, ... 
    maxCacheSize, stepType, gap_check, rand_seed, datasetName);


% % Note that in our experiments, we used an APT toolbox designed for our cluster
% % (https://github.com/iXce/APT). Here is the exact code that we ran:
%
% APT_params();
% global APT_PARAMS;
% APT_PARAMS.exec_name = 'horse_BCFW_hybrid';
% APT_PARAMS.host_name = { 'mem001' };
% APT_compile( 'horse_BCFW_hybrid.m' );
% 
% APT_run('horse_BCFW_hybrid.m', ...
%     {dataPath}, {resultPath}, lambda, gap_threshold, num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, maxCacheSize, stepType, gap_check, rand_seed, datasetName, ...
%     'Memory', 20*1024, 'CombineArgs', 1, 'ClusterID', 2 );
