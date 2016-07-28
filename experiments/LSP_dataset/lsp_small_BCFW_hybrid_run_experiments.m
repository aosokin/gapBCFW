%% setup paths
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

dataPath = fullfile(rootPath, 'data', 'LSP' );
resultPath = fullfile(rootPath, 'results' , 'lsp_small_BCFW_hybrid' );

%% run the jobs
lambda             = {10, 100, 1000};
gap_threshold      = {10^-8};
num_passes         = {50000};
time_budget        = {2880}; % 48h
sample             = {'uniform', 'gap'};
useCache           = {1, 0};
cacheNu            = {0.01};
cacheFactor        = {0.25};
maxCacheSize       = {20};
stepType           = {1,2};
gap_check          = {10};
rand_seed          = num2cell(6:10);


task_launcher('lsp_small_BCFW_hybrid', ...
    {dataPath},{resultPath}, lambda, gap_threshold, num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, maxCacheSize, stepType, gap_check, rand_seed);

% % Note that in our experiments, we used a toolbox designed for our cluster
% % named APT (https://github.com/iXce/APT), here is the exact code that run
% % for our cluster:
%
% APT_params();
% global APT_PARAMS;
% APT_PARAMS.exec_name = 'ocr_large_BCFW_hybrid';
% APT_PARAMS.host_name = { 'mem001' };
% APT_compile( 'lsp_small_BCFW_hybrid.m' );
% APT_run('lsp_small_BCFW_hybrid.m', ...
%     {dataPath},{resultPath}, lambda, gap_threshold, num_passes, time_budget, sample, useCache, cacheNu, cacheFactor, maxCacheSize, stepType, gap_check, rand_seed, ...
%     'Memory', 10*1024, 'CombineArgs', 1, 'ClusterID', 2 );
