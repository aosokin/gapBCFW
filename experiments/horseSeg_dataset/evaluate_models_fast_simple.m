function info = evaluate_models_fast_simple( progress, param, lambda, patterns_train, labels_train, patterns_test, labels_test )
%evaluate_models_fast_simple is a n even faster version of eveluate_models_fast.m
% evaluate_models_fast_simple does not compute loss on the training set so is better suited for large datasets

%% prepare output
info = struct;       

%% compute objective on the training set
    % reuse all gap computations done at the training time; augment them with first and last models
    info.train.gap = progress.verifiedGap;
    info.train.modelId = progress.verifiedGapModelId;
    if info.train.modelId(1) ~= 1
        fprintf('Computing the gap for the first model\n');
        info.train.gap = [nan; info.train.gap(:)];
        info.train.gap(1) = duality_gap( param, param.oracleFn, progress.models{1}, lambda );
        info.train.modelId = [1; info.train.modelId(:)];
    end
    if info.train.modelId(end) ~= numel(progress.models)
        fprintf('Computing the gap for the last model\n');
        info.train.gap = [info.train.gap(:); nan];
        info.train.gap(end) = duality_gap( param, param.oracleFn, progress.models{end}, lambda );
        info.train.modelId = [info.train.modelId(:); numel(progress.models)];
    end
    info.train.time = progress.time(info.train.modelId);
    info.verifiedGap = progress.verifiedGap;
    info.verifiedGapModelId = progress.verifiedGapModelId;
    info.verifiedGapNumPasses = progress.verifiedGapNumPasses;
    info.train.numPasses = progress.numPasses(info.train.modelId);

    num_models_evaluation = numel(info.train.modelId);                                        

    info.train.dual = nan(num_models_evaluation, 1);
    info.train.primal = nan(num_models_evaluation, 1);
    
    for i_model = 1 : num_models_evaluation
         model = progress.models{ info.train.modelId(i_model) };
         info.train.dual(i_model) = -objective_function(model.w, model.ell, lambda); % dual value -equation (4)
         info.train.primal(i_model) = info.train.dual(i_model) + info.train.gap(i_model);
    end
    
% %% compute error on the training set
% for i_model = 1 : num_models_evaluation
%     fprintf( 'Train error: model %d of %d\n', i_model, num_models_evaluation );
%     model = progress.models{ info.train.modelId(i_model) };
% 
%     % loss on train set
%     avg_loss = 0;
%     for i=1:numel(patterns_train)
%         ypredict = param.oracleFn(param, model, patterns_train{i}); % standard prediction as don't give label as input
%         avg_loss = avg_loss + param.lossFn(param, labels_train{i}, ypredict);
%     end
%     avg_loss = avg_loss / numel(patterns_train);
%     
%     info.train.error( i_model ) = avg_loss;
% end

%% compute error on the test set
subsampleStep = 1;
if numel(patterns_test) > numel(patterns_train)
    subsampleStep = floor(numel(patterns_test) / numel(patterns_train));
end

info.test.modelId = info.train.modelId( [1 : subsampleStep : num_models_evaluation-1, num_models_evaluation] ); % make sure that the last model is taken
info.test.time =  progress.time(info.test.modelId);
info.test.numPasses = progress.numPasses(info.test.modelId);
info.test.error = nan( numel(info.test.modelId), 1 );
for i_model =  1 : numel(info.test.modelId)
    fprintf( 'Test error: model %d of %d\n', i_model, numel(info.test.modelId) );

   model = progress.models{ info.test.modelId(i_model) };
    % loss on the test set
    avg_loss = 0;
    for i=1:numel(patterns_test)
        ypredict = param.oracleFn(param, model, patterns_test{i}); % standard prediction as don't give label as input
        avg_loss = avg_loss + param.lossFn(param, labels_test{i}, ypredict);
    end
    avg_loss = avg_loss / numel(patterns_test);
    
    info.test.error( i_model ) = avg_loss;
end

end

