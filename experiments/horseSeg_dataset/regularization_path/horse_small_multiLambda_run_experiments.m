%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

dataPath = fullfile(rootPath, 'data', 'horseSeg' );
resultPath = fullfile(rootPath, 'results' , 'horseSeg_small_regPath' );

%% prepare the job
dataPath_cell = {};
resultPath_cell = {};
lambda_grid = {};
gap_threshold = {};
num_passes = {};
time_budget = {};
sample = {};
useCache = {};
cacheNu = {};
cacheFactor = {};
maxCacheSize = {};
stepType = {};
gap_check = {};
warm_start_type = {};
datasetName = {};
rand_seed = {};
rand_seeds_all = 6 : 10;
time_limit = 10*60;

% job 1 : 
% - warm start
% NumPasses = GAP + pFW + Cache + primal warm start
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    lambda_grid{end+1} = 'factor2exp15to-15';
    gap_threshold{end+1} = 0.1;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 0;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 1;
    gap_check{end+1} = 10;
    warm_start_type{end+1} = 'keep_primal';
    rand_seed{end+1} = i_seed;
    datasetName{end+1} = 'horseSeg_small';
end

% job 2 : 
% - warm start
% Time =  GAP+BCFW+No Cache + primal warm start
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    lambda_grid{end+1} = 'factor2exp15to-15';
    gap_threshold{end+1} = 0.1;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 1;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 2;
    gap_check{end+1} = 10;
    warm_start_type{end+1} = 'keep_primal';
    rand_seed{end+1} = i_seed;
    datasetName{end+1} = 'horseSeg_small';
end

% job 3 : 
% - cold start
% NumPasses = GAP + pFW + Cache + no warm start
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    lambda_grid{end+1} = 'factor2exp15to-15';
    gap_threshold{end+1} = 0.1;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 0;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 1;
    gap_check{end+1} = 10;
    warm_start_type{end+1} = 'none';
    rand_seed{end+1} = i_seed;
    datasetName{end+1} = 'horseSeg_small';
end

% job 4 : 
% - cold start
% Time =  GAP+BCFW+No Cache + no warm start
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    lambda_grid{end+1} = 'factor2exp15to-15';
    gap_threshold{end+1} = 0.1;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 1;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 2;
    gap_check{end+1} = 10;
    warm_start_type{end+1} = 'none';
    rand_seed{end+1} = i_seed;
    datasetName{end+1} = 'horseSeg_small';
end

%% run the jobs
task_launcher('horse_multiLambda', ...
    dataPath_cell, resultPath_cell, lambda_grid, gap_threshold, ...
    num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, ...
    maxCacheSize, stepType, gap_check, warm_start_type, rand_seed, datasetName);


% % Note that in our experiments, we used an APT toolbox designed for our cluster
% % (https://github.com/iXce/APT). Here is the exact code that we ran:
%
% APT_params();
% global APT_PARAMS;
% APT_PARAMS.exec_name = 'horse_multiLambda';
% APT_PARAMS.host_name = { 'mem001' };
% APT_compile( 'horse_multiLambda.m' );
% 
% APT_run('horse_multiLambda.m', ...
%     dataPath_cell, resultPath_cell, lambda_grid, gap_threshold, ...
%     num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, ...
%     maxCacheSize, stepType, gap_check, warm_start_type, rand_seed, datasetName, ...
%     'Memory', 10*1024, 'CombineArgs', 0, 'ClusterID', 2);


