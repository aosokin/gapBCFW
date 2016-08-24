%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

dataPath = fullfile(rootPath, 'data', 'OCR' );
resultPath = fullfile(rootPath, 'results' , 'ocr_large_regPath' );

%% prepare the job
dataPath_cell = {};
resultPath_cell = {};
num_passes = {};
time_budget = {};
sample = {};
useCache = {};
cacheNu = {};
cacheFactor = {};
maxCacheSize = {};
stepType = {};
rand_seed = {};
AplusB = {};
A = {};
true_gap = {};
eps_reg_path = {};
gap_check = {};
datasetName = {};
rand_seeds_all = 6 : 10;
time_limit = 24*60;

% job 1 : 
% - exact gap
% NumPasses = GAP + pFW + Cache + (A+B=1) + (A=0.9)
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 1;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 2;
    rand_seed{end+1} = i_seed;
    AplusB{end+1} = 1;
    A{end+1} = 0.9;
    true_gap{end+1} = 1;
    eps_reg_path{end+1} = 0.1;
    gap_check{end+1} = 10;
    datasetName{end+1} = 'large';
end

% job 2 : 
% - exact gap
% Time = GAP + BCFW + No Cache +  (A+B=1)+ (A=0.9)
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 0;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 1;
    rand_seed{end+1} = i_seed;
    AplusB{end+1} = 1;
    A{end+1} = 0.9;
    true_gap{end+1} = 1;
    eps_reg_path{end+1} = 0.1;
    gap_check{end+1} = 10;
    datasetName{end+1} = 'large';
end

% job 3 : 
% - heuristic reg path 
% NumPasses = GAP + pFW + Cache + (A+B=1) (accurate enough) + (A=0.7)
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 1;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 2;
    rand_seed{end+1} = i_seed;
    AplusB{end+1} = 1;
    A{end+1} = 0.7;
    true_gap{end+1} = 0;
    eps_reg_path{end+1} = 0.1;
    gap_check{end+1} = 10000;
    datasetName{end+1} = 'large';
end

% job 4 : 
% - heuristic reg path 
% Time = GAP + BCFW + No Cache +  (A+B=1) (accurate enough) + (A=0.7)
for i_seed = rand_seeds_all
    dataPath_cell{end+1} = dataPath;
    resultPath_cell{end+1} = resultPath;
    num_passes{end+1} = 200000;
    time_budget{end+1} = time_limit;
    sample{end+1} = 'gap';
    useCache{end+1} = 0;
    cacheNu{end+1} = 0.01;
    cacheFactor{end+1} = 0.25;
    maxCacheSize{end+1} = 50;
    stepType{end+1} = 1;
    rand_seed{end+1} = i_seed;
    AplusB{end+1} = 1;
    A{end+1} = 0.7;
    true_gap{end+1} = 0;
    eps_reg_path{end+1} = 0.1;
    gap_check{end+1} = 10000;
    datasetName{end+1} = 'large';
end


%% run the jobs
task_launcher('ocr_regPath', ...
    dataPath_cell, resultPath_cell, time_budget, ...
    true_gap, A, AplusB, stepType, eps_reg_path, gap_check, sample, useCache, ...
    num_passes, maxCacheSize, cacheNu, cacheFactor, rand_seed, datasetName);


% % Note that in our experiments, we used an APT toolbox designed for our cluster
% % (https://github.com/iXce/APT). Here is the exact code that we ran:
%
% APT_params();
% global APT_PARAMS;
% APT_PARAMS.exec_name = 'ocr_regPath';
% APT_PARAMS.host_name = { 'mem001' };
% APT_compile( 'ocr_regPath.m' );
% 
% APT_run('ocr_regPath.m', ...
%     dataPath_cell, resultPath_cell, time_budget, ...
%     true_gap, A, AplusB, stepType, eps_reg_path, gap_check, sample, useCache, ...
%     num_passes, maxCacheSize, cacheNu, cacheFactor, rand_seed, datasetName, ...
%     'Memory', 18*1024, 'CombineArgs', 0, 'ClusterID', 2);

