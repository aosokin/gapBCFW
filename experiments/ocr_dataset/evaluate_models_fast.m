function info = evaluate_models_fast( progress, param, lambda, patterns_train, labels_train, patterns_test, labels_test )
%evaluate_models_fast is a faster version of eveluate_models.m that tries to save computations

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
    

% do it only for 150 models maximum
num_models = numel( progress.models );
step_size = ceil(num_models/150);
info.step_size = step_size;

for i_model = 1 : step_size : num_models
    fprintf( 'Train error: model %d of %d\n', i_model, num_models );
    model = progress.models{ i_model };

    % loss on the train set
    avg_loss = 0;
    for i=1:numel(patterns_train)
        ypredict = param.oracleFn(param, model, patterns_train{i}); % standard prediction as don't give label as input
        avg_loss = avg_loss + param.lossFn(param, labels_train{i}, ypredict);
    end
    avg_loss = avg_loss / numel(patterns_train);
    
    info.train.error( i_model ) = avg_loss;
end


%% compute error on the test set
for i_model = 1 : step_size :  num_models
    fprintf( 'Test error: model %d of %d\n', i_model, num_models );

    model = progress.models{ i_model };

    % loss on test set
    avg_loss = 0;
    for i=1:numel(patterns_test)
        ypredict = param.oracleFn(param, model, patterns_test{i}); % standard prediction as don't give label as input
        avg_loss = avg_loss + param.lossFn(param, labels_test{i}, ypredict);
    end
    avg_loss = avg_loss / numel(patterns_test);
    
    info.test.error( i_model ) = avg_loss;
end


end

