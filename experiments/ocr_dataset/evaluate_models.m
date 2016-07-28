function info = evaluate_models( models, param, lambda, patterns_train, labels_train, patterns_test, labels_test )
% evaluate_models computes primal/dual/gap, train error, test error for multiple models

%% prepare output
num_models = length(models);                                        
info = struct;
info.train.primal = nan( num_models, 1 );
info.train.dual = nan( num_models, 1 );
info.train.gap = nan( num_models, 1 );
info.train.error = nan( num_models, 1 );

info.test.error = nan( num_models, 1 );

%% compute objective on the training set
for i_model = 1 : num_models
    fprintf( 'Train objective: model %d of %d\n', i_model, num_models );

    model = models{ i_model };
    
    [gap, w_s, ell_s] = duality_gap( param, param.oracleFn, model, lambda );
    dual = -objective_function(model.w, model.ell, lambda); % dual value -equation (4)
    primal = dual+gap; % a cheaper alternative to get the primal value
    
    info.train.gap( i_model ) = gap;
    info.train.dual( i_model ) = dual;
    info.train.primal( i_model ) = primal;
end

%% compute error on the training set
for i_model = 1 : num_models
    fprintf( 'Train error: model %d of %d\n', i_model, num_models );
    model = models{ i_model };

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
subsampleStep = 1;
if numel(patterns_test) > numel(patterns_train)
    subsampleStep = floor(numel(patterns_test) / numel(patterns_train));
end
    
for i_model = 1 : subsampleStep : num_models
    fprintf( 'Test error: model %d of %d\n', i_model, num_models );

    model = models{ i_model };

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

