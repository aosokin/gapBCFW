function info = evaluate_regPath_models( param, progress )

num_lambda = numel(progress);
info = struct;

for i_lambda = 1 : num_lambda
    lambda = progress{i_lambda}.lambda;
    fprintf('Evaluating model for lambda=%f (%d of %d)\n', lambda, i_lambda, num_lambda);
    
    info(i_lambda).lambda = lambda;
    
    model = progress{i_lambda}.models{end};
    info(i_lambda).gap = duality_gap( param, param.oracleFn, model, lambda );
    
    info(i_lambda).time = progress{i_lambda}.time(end);
    info(i_lambda).numPasses = progress{i_lambda}.numPasses(end);
    info(i_lambda).numModels = numel( progress{i_lambda}.numPasses );
end


end

