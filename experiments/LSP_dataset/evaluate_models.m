function info = evaluate_models( models, param, lambda, patterns_train, labels_train, patterns_test, labels_test, progress )
                                            
%% prepare output
num_models = length(models);                                        
info = struct;
info.train.primal = nan( num_models, 1 );
info.train.dual = nan( num_models, 1 );
info.train.gap = nan( num_models, 1 );
info.train.error = nan( num_models, 1 );

info.test.error = nan( num_models, 1 );


%% compute objective on the training set (not necessary because we already check gaps)
if nargin==6
    for i_model = 1 : 10 :  num_models
        fprintf( 'Train objective: model %d of %d\n', i_model, num_models );

        model = models{ i_model };

        [gap, w_s, ell_s] = duality_gap( param, param.oracleFn, model, lambda );
        dual = -objective_function(model.w, model.ell, lambda); % dual value -equation (4)
        primal = dual+gap; % a cheaper alternative to get the primal value

        info.train.gap( i_model ) = gap;
        info.train.dual( i_model ) = dual;
        info.train.primal( i_model ) = primal;
    end
else
    for i_model = progress.verifiedGapModelId.'
        fprintf( 'Train objective: model %d of %d\n', i_model, num_models );

        model = models{ i_model };
        
        idx_gap = progress.verifiedGapModelId==i_model;
        gap = progress.verifiedGap(idx_gap);
        dual = -objective_function(model.w, model.ell, lambda); % dual value -equation (4)
        primal = dual+gap; % a cheaper alternative to get the primal value

        info.train.gap( i_model ) = gap;
        info.train.dual( i_model ) = dual;
        info.train.primal( i_model ) = primal;
    end
end

%% compute error on the training set

% do it only for 10 models maximum

step_size = ceil(num_models/3);
info.step_size = step_size;

% for i_model = 1 : step_size : num_models
%     fprintf( 'Train error: model %d of %d\n', i_model, num_models );
%     model = models{ i_model };
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
%     
%    % info.train.perf_paper(i_model) = test_model_comparisonpaper(model, param, patterns_train);
% end


%% compute error on the test set
for i_model = 1 : step_size :  num_models
    fprintf( 'Test error: model %d of %d\n', i_model, num_models );

    model = models{ i_model };

    % loss on train set
    avg_loss = 0;
    for i=1:5:numel(patterns_test)
        ypredict = param.oracleFn(param, model, patterns_test{i}); % standard prediction as don't give label as input
        avg_loss = avg_loss + param.lossFn(param, labels_test{i}, ypredict);
    end
    avg_loss = avg_loss / numel(patterns_test);
    
    info.test.error( i_model ) = avg_loss;
    
    % now computing the error as reported in the paper of Chen and Yuille
    
    % getting performance from the paper
    %info.test.perf_paper(i_model) = test_model_comparisonpaper(model, param, patterns_test);
        
end

end

