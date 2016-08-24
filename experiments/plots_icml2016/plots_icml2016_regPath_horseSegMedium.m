% Figure 8c

%% Get the default result path
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
    rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

resultPath = fullfile(rootPath, 'results' , 'horseSeg_medium_regPath' );

%% general plot settings
saveFig     = false; %save to true if you are happy with the figures
pathSaveFig = './figures';
dataSet     = 'horseMedium';
revertAxis   = true;

firstRun   = false; % if you don't know scale of things in advance

lambdaTarget = [2.^(-15:0.5:-5), 2.^(-5:30)]; 
xlimlambda  = [10^-3 10^6];
maxTenPower = ceil(log(max(xlimlambda))/log(10));
minTenPower = floor(log(min(xlimlambda))/log(10));
xtickLambda = 10.^(minTenPower:2:maxTenPower);

ylimtime    = [10^-2 30];
maxTenPower = ceil(log(max(ylimtime))/log(10));
minTenPower = floor(log(min(ylimtime))/log(10));
ytickTime = 10.^(minTenPower:1:maxTenPower);

ylimnumpass = [10^0 10^3.8];
maxTenPower = ceil(log(max(ylimnumpass))/log(10));
minTenPower = floor(log(min(ylimnumpass))/log(10));
ytickNumPasses = 10.^(minTenPower:1:maxTenPower);

timeFormat = 'h';

idFig=1;
rng(1);
figure_size = [450, 300];


%% prepare file names for plotting
resultFiles_bcfw = cell(4, 1);
resultFiles_bcPfw = cell(4, 1);

% exact path
setupTemplateName = fullfile(resultPath, ['results_regPath', ...
    '_numPasses','%s',...
    '_seed','%s',...
    '_gapCheck','%s',...
    '_sample_','%s',...
    '_useCache','%s',...
    '_cacheNu','%s',...
    '_cacheFactor','%s',...
    '_stepType','%s',...
    '_timeBudget','%s',...
    '_maxCacheSize','%s',...
    '_A','%s',...
    '_AplusB','%s',...
    '_epsRegPath','%s',...
    '_trueGap','%s',...
    '.mat' ] );
gap_check_value = 10;
eps_reg_path_value = 0.1;
AplusB_value = 1;
A_value = 0.9;
true_gap_value = 1;
time_budget_value = 24*60;
num_passes_value = 200000;
maxCacheSize_value = 50;
cacheNu_value = 0.01;
cacheFactor_value = 0.25;
sample_value = 'gap';

curFileTemplate = sprintf( setupTemplateName, num2str(num_passes_value), '%s', num2str(gap_check_value), ...
    sample_value, '%s', num2str(cacheNu_value), num2str(cacheFactor_value), '%s', ...
    num2str(time_budget_value), num2str(maxCacheSize_value), num2str(A_value), num2str(AplusB_value), ...
    num2str(eps_reg_path_value), num2str(true_gap_value) );

useCache_value = 0;
stepType_value = 1;
resultFiles_bcfw{2} = sprintf( curFileTemplate, '%s', num2str(useCache_value), num2str(stepType_value) );

useCache_value = 1;
stepType_value = 2;
resultFiles_bcPfw{2} = sprintf( curFileTemplate, '%s', num2str(useCache_value), num2str(stepType_value) );

% heuristic path
true_gap_value = 0;
A_value = 0.7;
gap_check_value = 10000;

curFileTemplate = sprintf( setupTemplateName, num2str(num_passes_value), '%s', num2str(gap_check_value), ...
    sample_value, '%s', num2str(cacheNu_value), num2str(cacheFactor_value), '%s', ...
    num2str(time_budget_value), num2str(maxCacheSize_value), num2str(A_value), num2str(AplusB_value), ...
    num2str(eps_reg_path_value), num2str(true_gap_value) );

useCache_value = 0;
stepType_value = 1;
resultFiles_bcfw{1} = sprintf( curFileTemplate, '%s', num2str(useCache_value), num2str(stepType_value) );

useCache_value = 1;
stepType_value = 2;
resultFiles_bcPfw{1} = sprintf( curFileTemplate, '%s', num2str(useCache_value), num2str(stepType_value) );

% grid search: no warm start
setupTemplateName = fullfile(resultPath, ['results_gridSearch', ...
    '_lambdaGrid_', '%s', ...
    '_gapThreshold',  '%s', ...
    '_warmStart',  '%s', ...
    '_numPasses',  '%s', ...
    '_timeBudget', '%s', ...
    '_sample_',  '%s', ...
    '_useCache',  '%s', ...
    '_cacheNu',  '%s', ...
    '_cacheFactor',  '%s', ...
    '_maxCacheSize',  '%s', ...
    '_stepType', '%s', ...
    '_gapCheck',  '%s', ...
    '_seed',  '%s', ...
    '.mat' ]);
lambda_value = 'factor2exp15to-15';
gap_threshold_value = 0.1;
num_passes_value = 200000;
time_budget_value = 24 * 60;
sample_value = 'gap';
cacheNu_value = 0.01;
cacheFactor_value = 0.25;
maxCacheSize_value = 50;
gap_check_value = 10;
warm_start_type_value = 'none';

curFileTemplate = sprintf( setupTemplateName, lambda_value, num2str(gap_threshold_value), warm_start_type_value, num2str(num_passes_value), ...
    num2str(time_budget_value), sample_value, '%s', num2str(cacheNu_value), num2str(cacheFactor_value), ...
    num2str(maxCacheSize_value), '%s', num2str(gap_check_value), '%s' );

useCache_value = 0;
stepType_value = 1;
resultFiles_bcfw{3} = sprintf( curFileTemplate, num2str(useCache_value), num2str(stepType_value), '%s' );

useCache_value = 1;
stepType_value = 2;
resultFiles_bcPfw{3} = sprintf( curFileTemplate, num2str(useCache_value), num2str(stepType_value), '%s' );

% grid search: no warm start
warm_start_type_value = 'keep_primal';
curFileTemplate = sprintf( setupTemplateName, lambda_value, num2str(gap_threshold_value), warm_start_type_value, num2str(num_passes_value), ...
    num2str(time_budget_value), sample_value, '%s', num2str(cacheNu_value), num2str(cacheFactor_value), ...
    num2str(maxCacheSize_value), '%s', num2str(gap_check_value), '%s' );

useCache_value = 0;
stepType_value = 1;
resultFiles_bcfw{4} = sprintf( curFileTemplate, num2str(useCache_value), num2str(stepType_value), '%s' );

useCache_value = 1;
stepType_value = 2;
resultFiles_bcPfw{4} = sprintf( curFileTemplate, num2str(useCache_value), num2str(stepType_value), '%s' );


%% load data
seed_values = 6 : 10;
data_bcfw = cell(4, 1);
data_bcPfw = cell(4, 1);
for i_seed = 1 : length(seed_values)
    for i_method = 1 : max( length(resultFiles_bcfw), length(resultFiles_bcPfw) )
        data_bcfw{i_method}{i_seed} = load( sprintf(resultFiles_bcfw{i_method}, num2str(seed_values(i_seed)) ), ...
            'info');
        
        data_bcPfw{i_method}{i_seed} = load( sprintf(resultFiles_bcPfw{i_method}, num2str(seed_values(i_seed)) ), ...
            'info');
    end
end

%% plot everything
for plotNumber=1:4
    
    switch plotNumber
        case 1
            criterion = 'time';
            methodName = 'bcfwGap';
        case 2
            criterion = 'time';
            methodName = 'bcPfwGapCache';
        case 3
            criterion = 'numpass';
            methodName = 'bcfwGap';
        case 4
            criterion = 'numpass';
            methodName = 'bcPfwGapCache';
    end
    
    nameSave = [dataSet,'_',methodName,'_', criterion];
    
    plotLine      = {'Heur', 'Exact', 'GridWarm', 'GridCold'}; % here we want exact but could be Grid or Heur
    
    result = [];
    for p=1:numel(plotLine)
        tmp_res = struct();
        switch methodName
            case 'bcfwGap'
                res = data_bcfw{p};
            case 'bcPfwGapCache'
                res = data_bcPfw{p};
        end
        
        tmp_res.xValues = cell(0);
        tmp_res.yValues = cell(0);
        for i_seed = 1 : length(seed_values)
            tmp_res.xValues{i_seed} = cat(1, res{i_seed}.info.lambda);
            
            switch criterion
                case 'time'
                    tmp_res.yValues{i_seed} = cumsum(cat(1, res{i_seed}.info.time));
                    switch timeFormat
                        case 'h'
                            tmp_res.yValues{i_seed} = tmp_res.yValues{i_seed} ./ 3600;
                        case 'm'
                            tmp_res.yValues{i_seed} = tmp_res.yValues{i_seed} ./ 60;
                        case 's'
                            tmp_res.yValues{i_seed} = tmp_res.yValues{i_seed};
                    end
                    
                case 'numpass'
                    tmp_res.yValues{i_seed} = cumsum(cat(1, res{i_seed}.info.numPasses));
            end
        end
        
        % filter out results that you don't want (here based on A
        % parameters)
        if strcmp(plotLine{p},'Exact')
            tmp_res.name = '\epsilon-approximate Path';
            tmp_res.color = [41,208,208]./255;
            tmp_res.marker = {'Marker', 'square', 'MarkerEdgeColor', 'k'};
            tmp_res.linestyle = {'LineStyle', '-', 'LineWidth', 1.5};
        end
        if strcmp(plotLine{p},'Heur')
            tmp_res.name = 'Heuristic Path';
            tmp_res.color = [41,208,208]./255;
            tmp_res.marker = {'Marker', 'o', 'MarkerEdgeColor', 'k'};
            tmp_res.linestyle = {'LineStyle', '-', 'LineWidth', 1.5};
        end
        
        if strcmp(plotLine{p}, 'GridWarm')
            tmp_res.name = 'Grid: warm start';
            tmp_res.color = [255,0,0]./255;
            tmp_res.marker = {'Marker', 'o', 'MarkerEdgeColor', 'k'};
            tmp_res.linestyle = {'LineStyle', '-', 'LineWidth', 1.5};
            
        end
        if strcmp(plotLine{p}, 'GridCold')
            tmp_res.name = 'Grid: no warm start';
            tmp_res.color = [255,0,0]./255;
            tmp_res.marker = {'Marker', 'square', 'MarkerEdgeColor', 'k'};
            tmp_res.linestyle = {'LineStyle', '-', 'LineWidth', 1.5};
        end
        
        if isempty(result)
            result = tmp_res;
        else
            result = [result tmp_res];
        end
        
    end
    
    
    
    h = figure(idFig);
    clf;
    hold on;
    idFig = idFig+1;
    hLeg = plotMultipleMethodsWithErrorsBars(result, lambdaTarget);
    xlabel('Effective passes over data');
    ylabel('Expected improvement');
    cur_pos = get(gcf,'Position');
    set(h, 'Position', [cur_pos(1), cur_pos(2), figure_size(1), figure_size(2)]);
    hold off
    
    
    
    % change limits of the plot
    switch criterion
        case 'time'
            set(gca, 'YScale', 'log');
            set(gca, 'XScale', 'log');
            xlabel('Value of \lambda');
            ylabel(sprintf('Time in (%s)', timeFormat));
            set(gca, 'FontSize', 12);
            if ~firstRun
                xlim(xlimlambda);          % directly change here
                ylim(ylimtime);
                set(gca,'XTick', xtickLambda);
                set(gca,'YTick', ytickTime);
            end
            
            if revertAxis
                set(gca,'xdir','reverse');
                set(hLeg,'Location','NorthWest','box','off');
            else
                set(hLeg,'Location','SouthWest','box','off');
            end
            
        case 'numpass'
            set(gca, 'YScale', 'log');
            set(gca, 'XScale', 'log');
            xlabel('Value of \lambda');
            ylabel('Num passes');
            set(gca, 'FontSize', 12);
            if ~firstRun
                xlim(xlimlambda);          % directly change here
                ylim(ylimnumpass);
                set(gca,'XTick', xtickLambda);
                set(gca,'YTick', ytickNumPasses);
            end
            
            if revertAxis
                set(gca,'xdir','reverse');
                set(hLeg,'Location','NorthWest','box','off');
            else
                set(hLeg,'Location','SouthWest','box','off');
            end
    end
    
    if saveFig
        title('');
        % remove legend?
        name_save = fullfile(pathSaveFig, nameSave);
        name_crop = fullfile(sprintf('%s-crop.pdf', name_save));
        hgexport(h, name_save, hgexport('factorystyle'), 'Format', 'pdf');
        system(sprintf('pdfcrop %s.pdf', name_save));
        system(sprintf('mv %s %s.pdf', name_crop, name_save));
    end
    
end
