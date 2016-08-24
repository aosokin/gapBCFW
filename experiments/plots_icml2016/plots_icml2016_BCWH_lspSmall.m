% Figure 6 LSP-small

%% init
methodDisplayName = cell(0,0);
methodFile = cell(0,0);
methodLine = cell(0,0);
methodColor = cell(0,0);
methodMarker = cell(0,0);

formatResult = ['results_BCFW_hybrid', ...
    '_lambda',     '%s', ...
    '_gapThreshold','%s', ...
    '_numPasses','%s',...
    '_timeBudget', '%s' ...
    '_sample_', '%s', ...
    '_useCache', '%s', ...
    '_cacheNu', '%s', ...
    '_cacheFactor', '%s', ...
    '_maxCacheSize', '%s', ...
    '_stepType', '%s', ...
    '_gapCheck', '%s', ...
    '_seed','%s.mat'];

%% Get the default result path
rootPath = fileparts( mfilename('fullpath') );
while ~exist( fullfile( rootPath, 'setup_BCFW.m' ), 'file' )
     rootPath = fileparts( rootPath );
end
run( fullfile( rootPath, 'setup_BCFW' ) );

resultPath = fullfile(rootPath, 'results' , 'lsp_small_BCFW_hybrid' );

%% parameters for cache
lambdaValues       = {10,100,1000};
cachingValues      = {0,1};
gapThresholdValues = {10^-8};
numPassesValues    = {50000};
timeBudgetValues   = {2880}; % 48h
samplingValues     = {'uniform', 'gap'};
cacheNus           = {0.01};
cacheFactors       = {0.25};
maxCacheSizeValues = {20};
stepTypeValues     = {1,2};
gapCheckValues     = {10};
seedValues         = num2cell(6:10);

%% parameters for plot
rng(1);

dataSet           = 'lspSmall'; % don't forget to also change the result path and everything that depend on the dataset
pathSaveFig       = './figures';
saveFig           = false; % flag to save figure
iFig              = 1; % index of the first figure
alpha             = 0; % quantile parameter
firstRun          = false; % if this is on then the xlim and ylim are ignored (maybe useful if you don't know where should be the limits)
displayTrainError = false; % flag for train error
displayTestError  = false;

% limit size parameters (note that you can specify a table if you want
% different scale for different lambda)

% first plot Gap vs Passes
nPtsMarkersGapVsPass = 3;
xPassTargetVsGap     = 0:35:700;
xLimGapVsPass        = [0, 700];
yLimGapVsPass        = [1.5*10^-2, 2;...
                        10^-3, 0.5;...
                        10^-4, 10^-1];

% second plot Gap vs Time
nPtsMarkersTimeVsGap = 3;
timeFormat           = 'h';
xTimeTargetVsGap     = [0,0.1,0.3,0.5:1:49]; % this should have the same unit as time format
xLimGapVsTime        = [0, 48];
yLimGapVsTime        = yLimGapVsPass;

% third plot Test error vs Passes
nPtsMarkersTestErrorVsPass = 3;
xPassTargetVsTestError     = 1:20:700;
xLimTestErrorVsPass        = [0, 700];
yLimTestErrorVsPass        = [0.25, 0.6];

% Fourth plot (optional) Train error vs Passes
nPtsMarkersTrainErrorVsPass = 3;
xPassTargetVsTrainError     = 1:10:1000;
xLimTrainErrorVsPass        = [0, 500];
yLimTrainErrorVsPass        = [10^-4, 10^-5];

if saveFig && ~exist(pathSaveFig, 'dir')
    mkdir(pathSaveFig);
end

%% add the methods without cache
for tBV = 1:numel(timeBudgetValues)
    timeBudgetValue = timeBudgetValues{tBV};
    for gCV = 1:numel(gapCheckValues)
        gapCheckValue = gapCheckValues{gCV};
        for mCV=1:numel(maxCacheSizeValues)
            maxCacheSize = maxCacheSizeValues{mCV};
            for cF=1:numel(cacheFactors)
                cacheFactor = cacheFactors{cF};
                for cN=1:numel(cacheNus)
                    cacheNu = cacheNus{cN};
                    for nP=1:numel(numPassesValues)
                        numPasseValue = numPassesValues{nP};
                        for gTV=1:numel(gapThresholdValues)
                            gapThresholdValue = gapThresholdValues{gTV};
                            for sV=1:numel(samplingValues)
                                sampling = samplingValues{sV};
                                for m=1:numel(stepTypeValues)
                                    method = stepTypeValues{m};
                                    for c=1:numel(cachingValues)
                                        caching = cachingValues{c};
                                        [methodColor{end+1}, methodMarker{end+1}, methodLine{end+1}]  = getProperties(method, sampling, caching);
                                        methodDisplayName{end+1} = getName(method, sampling, caching);
                                        methodFile{end+1} = fullfile(resultPath, sprintf(formatResult, '%s', num2str(gapThresholdValue),...
                                            num2str(numPasseValue), num2str(timeBudgetValue), sampling, num2str(caching), num2str(cacheNu), ...
                                            num2str(cacheFactor), num2str(maxCacheSize), num2str(method), num2str(gapCheckValue), '%s'));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% do the plots
for l=1:numel(lambdaValues)
    %% LOOP OVER VALUES OF LAMBDA
    lambda = lambdaValues{l};
    %% plot gap vs oracle calls
    h = figure(iFig);
    iFig = iFig+1;
    clf;
    hold on
    numMethods    = length(methodDisplayName);
    numPassTarget = xPassTargetVsGap;
    nPtsMarkers   = nPtsMarkersGapVsPass;
    % fill the structure for results
    resultsGapVsOracle = struct();
    tablePlots  = [];
    for iMethod = 1 : numMethods
        if exist(sprintf(sprintf(methodFile{iMethod},num2str(lambda),num2str(6))), 'file')

            resultsGapVsOracle(iMethod).name      = methodDisplayName{iMethod};
            resultsGapVsOracle(iMethod).color     = methodColor{iMethod};
            resultsGapVsOracle(iMethod).linestyle = methodLine{iMethod};
            resultsGapVsOracle(iMethod).marker    = methodMarker{iMethod};

            % fill values to display
            curResults = cell(numel(seedValues),1);
            numPasses  = cell(numel(seedValues),1);
            for s=1:numel(seedValues)
                curRes = load( fullfile(sprintf(methodFile{iMethod}, num2str(lambda), num2str(seedValues{s}))));
                if isfield(curRes.train, 'numPasses')
                    numPasses{s} = curRes.train.numPasses;
                elseif isfield(curRes, 'numPasses')
                    numPasses{s} = curRes.numPasses;
                else
                    numPasses{s} = 1 : length(curResults.train.gap);
                end
                curResults{s} = curRes.train.gap;
            end

            resultsGapVsOracle(iMethod).xValues = numPasses;
            resultsGapVsOracle(iMethod).yValues = curResults;

        else
            warning(['Result file not found: ', methodFile{iMethod}]);

            resultsGapVsOracle(iMethod).name      = methodDisplayName{iMethod};
            resultsGapVsOracle(iMethod).color     = methodColor{iMethod};
            resultsGapVsOracle(iMethod).linestyle = methodLine{iMethod};
            resultsGapVsOracle(iMethod).marker    = methodMarker{iMethod};
            resultsGapVsOracle(iMethod).xValues = 1;
            resultsGapVsOracle(iMethod).yValues = 1;
        end
    end

    hLeg = plotMultipleMethodsWithErrorsBars(resultsGapVsOracle, numPassTarget, nPtsMarkers, alpha);

    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'linear');
    xlabel('Effective passes over data');
    ylabel('True duality gap');
    title( sprintf('%s, lambda = %0.0e', dataSet, lambda) );
    set(gca, 'FontSize', 9);
    set(hLeg, 'visible', 'off');

    if ~firstRun
        % check size of xlim
        if size(xLimGapVsPass,1)==1
            xlim(xLimGapVsPass);
        else
            xlim(xLimGapVsPass(l,:));
        end
        % check size of ylim
        if size(yLimGapVsPass,1)==1
            ylim(yLimGapVsPass);
        else
            ylim(yLimGapVsPass(l,:));
        end
    end
    cur_pos = get(h,'Position');
    set(h, 'Position', [cur_pos(1), cur_pos(2), 300, 200]);
%    grid on
    hold off
    if saveFig
        % not that you might first want to adapt the axis a bit
        title('');
        name_crop = fullfile(pathSaveFig, sprintf('gapVsPass_%s_lambda_%0.0e-crop.pdf', dataSet, lambda));
        name_save = fullfile(pathSaveFig, sprintf('gapVsPass_%s_lambda_%0.0e', dataSet, lambda));
        hgexport(h, name_save, hgexport('factorystyle'), 'Format', 'pdf');
        system(sprintf('pdfcrop %s.pdf', name_save));
        system(sprintf('mv %s %s.pdf', name_crop, name_save));
    end

    %% plot time vs dual gap

    h = figure(iFig);
    iFig = iFig+1;
    clf;
    hold on
    numMethods    = length(methodDisplayName);
    numTimeTarget = xTimeTargetVsGap;
    nPtsMarkers   = nPtsMarkersTimeVsGap;

    resultsGapVsTime = struct();
    for iMethod = 1 : numMethods
        if exist(sprintf(sprintf(methodFile{iMethod},num2str(lambda),num2str(seedValues{1}))), 'file')

            resultsGapVsTime(iMethod).name      = methodDisplayName{iMethod};
            resultsGapVsTime(iMethod).color     = methodColor{iMethod};
            resultsGapVsTime(iMethod).linestyle = methodLine{iMethod};
            resultsGapVsTime(iMethod).marker    = methodMarker{iMethod};

            % fill values to display
            curResults = cell(numel(seedValues),1);
            times      = cell(numel(seedValues),1);


            for s=1:numel(seedValues)
                curRes = load( fullfile(sprintf(methodFile{iMethod},num2str(lambda), num2str(seedValues{s}))));
                if isfield(curRes.train, 'time')
                    times{s} = transform_time(curRes.train.time, timeFormat);
                elseif isfield(curRes, 'time')
                    times{s} = transform_time(curRes.time, timeFormat);
                else
                    error('Time is not present! \n');
                end
                curResults{s} = curRes.train.gap;
            end

            resultsGapVsTime(iMethod).xValues = times;
            resultsGapVsTime(iMethod).yValues = curResults;

        else
            warning(['Result file not found: ', methodFile{iMethod}]);

            resultsGapVsTime(iMethod).name      = methodDisplayName{iMethod};
            resultsGapVsTime(iMethod).color     = methodColor{iMethod};
            resultsGapVsTime(iMethod).linestyle = methodLine{iMethod};
            resultsGapVsTime(iMethod).marker    = methodMarker{iMethod};
            resultsGapVsTime(iMethod).xValues   = 1;
            resultsGapVsTime(iMethod).yValues   = 1;
        end
    end

    hLeg = plotMultipleMethodsWithErrorsBars(resultsGapVsTime, numTimeTarget, nPtsMarkers, alpha);

    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'linear');
    xlabel(sprintf('Time (in %s)', timeFormat));
    ylabel('True duality gap');
    title( sprintf('%s, lambda = %0.0e', dataSet, lambda) );
    set(gca, 'FontSize', 9);
    set(hLeg, 'visible', 'off');

    if ~firstRun
        % check size of xlim
        if size(xLimGapVsTime,1)==1
            xlim(xLimGapVsTime);
        else
            xlim(xLimGapVsTime(l,:));
        end

        if size(yLimGapVsTime,1)==1
            ylim(yLimGapVsTime);
        else
            ylim(yLimGapVsTime(l,:));
        end
    end

    cur_pos = get(h, 'Position');
    set(h, 'Position', [cur_pos(1), cur_pos(2), 300, 200]);
    %set(h, 'Position', [cur_pos(1), cur_pos(2), 300, 200]);
    hold off
    if saveFig
        % not that you might first want to adapt the axis a bit
        title('');
        name_crop = fullfile(pathSaveFig, sprintf('gapVsTime_%s_lambda_%0.0e-crop.pdf', dataSet, lambda));
        name_save = fullfile(pathSaveFig, sprintf('gapVsTime_%s_lambda_%0.0e', dataSet, lambda));
        hgexport(h, name_save, hgexport('factorystyle'), 'Format', 'pdf');
        system(sprintf('pdfcrop %s.pdf', name_save));
        system(sprintf('mv %s %s.pdf', name_crop, name_save));
    end

    %% plot test error vs oracle calls

    if displayTestError
        h = figure(iFig);
        iFig = iFig+1;
        clf;
        hold on

        numPassTarget = xPassTargetVsTestError;
        numMethods    = length(methodDisplayName);
        nPtsMarkers   = nPtsMarkersTestErrorVsPass;

        resultsGapVsTestError = struct();

        for iMethod = 1 : numMethods
            curResults = cell(0,1);
            numPasses  = cell(0,1);
            curFile = sprintf(methodFile{iMethod}, num2str(lambda), num2str(seedValues{1}));
            if exist(curFile, 'file')
                resultsGapVsTestError(iMethod).name      = methodDisplayName{iMethod};
                resultsGapVsTestError(iMethod).color     = methodColor{iMethod};
                resultsGapVsTestError(iMethod).linestyle = methodLine{iMethod};
                resultsGapVsTestError(iMethod).marker    = methodMarker{iMethod};

                % fill values to display
                curResults = cell(numel(seedValues),1);
                numPasses  = cell(numel(seedValues),1);


                for s=1:numel(seedValues)
                    curRes = load( fullfile(sprintf(methodFile{iMethod}, num2str(lambda), num2str(seedValues{s}))));
                    if isfield(curRes.test, 'numPasses')
                        numPasses{s} = curRes.test.numPasses;
                    elseif isfield(curRes, 'train')
                        numPasses{s} = curRes.train.numPasses;
                    else
                        error('Result for numpasses is missing.');
                    end
                    curResults{s} = curRes.test.error;
                    fprintf('Error is %f for lambda %f and seed %f\n',curResults{s}(end), lambda, seedValues{s});


                end

                resultsGapVsTestError(iMethod).xValues = numPasses;
                resultsGapVsTestError(iMethod).yValues = curResults;
            else
                warning(['Result file not found: ', curFile]);
            end

        end

        hLeg = plotMultipleMethodsWithErrorsBars(resultsGapVsTestError, numPassTarget, nPtsMarkers, alpha);

        set(gca, 'YScale', 'log');
        set(gca, 'XScale', 'linear');
        xlabel('Effective passes over dataset');
        ylabel('Test error');
        title( sprintf('%s, lambda = %0.0e', dataSet, lambda) );
        set(hLeg, 'visible', 'off');
        set(gca, 'FontSize', 9);

        cur_pos = get(h,'Position');
        set(h, 'Position', [cur_pos(1), cur_pos(2), 300, 200]);


        if ~firstRun
            if size(xLimTestErrorVsPass,1)==1
                xlim(xLimTestErrorVsPass)
            else
                xlim(xLimTestErrorVsPass(l,:));
            end

            if size(yLimTestErrorVsPass,1)==1
                ylim(yLimTestErrorVsPass)
            else
                ylim(yLimTestErrorVsPass(l,:));
            end
        end

        hold off

        if saveFig
            % not that you might first want to adapt the axis a bit
            title('');
            name_save = fullfile(pathSaveFig, sprintf('testErrorVsPass_%s_lambda_%0.0e', dataSet, lambda));
            hgexport(h, name_save, hgexport('factorystyle'), 'Format', 'pdf');
            %saveas(h,name_save,'pdf');
            system(sprintf('pdfcrop %s.pdf', name_save));
            system(sprintf('mv %s %s.pdf', name_crop, name_save));
        end
    end

    if displayTrainError
        %% plot test error vs oracle calls
        h = figure(iFig);
        iFig = iFig+1;
        clf;
        hold on

        numPassTarget = xPassTargetVsTrainError;
        numMethods    = length(methodDisplayName);
        nPtsMarkers   = nPtsMarkersTrainErrorVsPass;

        resultsGapVsTrainError = struct();

        for iMethod = 1 : numMethods
            curResults = cell(0,1);
            numPasses  = cell(0,1);
            curFile = sprintf(methodFile{iMethod}, num2str(lambda), num2str(seedValues{1}));
            if exist(curFile, 'file')
                resultsGapVsTrainError(iMethod).name      = methodDisplayName{iMethod};
                resultsGapVsTrainError(iMethod).color     = methodColor{iMethod};
                resultsGapVsTrainError(iMethod).linestyle = methodLine{iMethod};
                resultsGapVsTrainError(iMethod).marker    = methodMarker{iMethod};

                % fill values to display
                curResults = cell(numel(seedValues),1);
                numPasses  = cell(numel(seedValues),1);

                for s=1:numel(seedValues)
                    curRes = load( fullfile(sprintf(methodFile{iMethod}, num2str(lambda), num2str(seedValues{s}))));
                    if isfield(curRes.train, 'numPasses')
                        numPasses{s} = curRes.train.numPasses;
                    elseif isfield(curRes, 'time')
                        numPasses{s} = curRes.numPasses;
                    else
                        error('Result for numpasses is missing.');
                    end
                    curResults{s} = curRes.train.error;
                end

                resultsGapVsTrainError(iMethod).xValues = numPasses;
                resultsGapVsTrainError(iMethod).yValues = curResults;
            else
                warning(['Result file not found: ', curFile]);
            end

        end

        hLeg = plotMultipleMethodsWithErrorsBars(resultsGapVsTrainError, numPassTarget, nPtsMarkers, alpha);

        set(gca, 'YScale', 'log');
        set(gca, 'XScale', 'linear');
        xlabel('Effective passes over dataset');
        ylabel('Train error');
        title( sprintf('%s, lambda = %0.0e', dataSet, lambda) );
        set(hLeg, 'visible', 'off');
        set(gca, 'FontSize', 9);

        cur_pos = get(h, 'Position');
        set(h, 'Position', [cur_pos(1), cur_pos(2), 300, 200]);


        if ~firstRun
            if size(xLimTrainErrorVsPass,1)==1
                xlim(xLimTrainErrorVsPass)
            else
                xlim(xLimTrainErrorVsPass(l,:));
            end

            if size(yLimTrainErrorVsPass,1)==1
                ylim(yLimTrainErrorVsPass)
            else
                ylim(yLimTrainErrorVsPass(l,:));
            end
        end

        hold off

        if saveFig
            % not that you might first want to adapt the axis a bit
            title('');
            name_save = fullfile(pathSaveFig, sprintf('trainErrorVsPass_%s_lambda_%0.0e', dataSet, lambda));
            hgexport(h, name_save, hgexport('factorystyle'), 'Format', 'pdf');
            system(sprintf('pdfcrop %s.pdf', name_save));
            system(sprintf('mv %s %s.pdf', name_crop, name_save));
        end

    end



end
