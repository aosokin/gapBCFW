function [gap, w_s, ell_s] = duality_gap( param, maxOracle, model, lambda )
% [gap, w_s, ell_s] = duality_gap( param, maxOracle, model, lambda )
% 
% Return the SVM duality gap for the implicit primal-dual pair given by
%   model.w and model.ell (w = A*\alpha; ell = b'*\alpha -- alpha is
%   implicit). See "Duality Gap" in Section 4 of ICML 2013 paper.
%
% This function is expensive, as it requires a full decoding pass over all
% examples (so it costs as much as n BCFW iterations). If the duality gap
% is checked regularly as a stopping criterion, then one can also use the
% returned w_s & ell_s quantities to make a batch Frank-Wolfe step and not
% waste this computation (see the update in Alg. 2 in the paper).
%   
% duality gap = lambda*(w-w_s)'*w - ell + ell_s
%
% ell_s = 1/n \sum_i ell(y_i, ystar_i) -- the average loss for the
%          *loss-augmented* predictions ystar_i
%
% w_s = 1/(lambda*n) \psi_i(ystar_i)

    loss = param.lossFn;
    phi = param.featureFn;
    
    w = model.w;
    ell = model.ell;
    
    n = numel(param.patterns);
    ystars = cell(n,1);
    for i=1:n
        % solve the loss-augmented inference for point i
        ystars{i} = maxOracle(param, model, param.patterns{i}, param.labels{i});
    end
    
    w_s = zeros(size(w));
    ell_s = 0; % note that we will normalize at the end to save rounding...
    for i=1:n
        psi_i = phi(param, param.patterns{i}, param.labels{i})-phi(param, param.patterns{i}, ystars{i});
        if isfield(param, 'cache_type_single') && param.cache_type_single
            % if cache is going to be single we need to truncate to single precision to avoid future in
            psi_i = single(psi_i);
            % need to cast back to double to have all the computations in double
            psi_i = double(psi_i);
        end
        w_s = w_s + psi_i;
        ell_s = ell_s + loss(param, param.labels{i}, ystars{i});
    end
      
    w_s = w_s/(lambda*n);
    ell_s = ell_s/n;
    
    
    % computing duality gap:
    % note that this expression is also the numerator of the line-search quotient solution
    gap = lambda* w'*(w - w_s) - ell + ell_s;
    
end % duality_gap
