function label = chain_conll_oracle(param, model, xi, yi)
% Do loss-augmented decoding on a given example (xi,yi) using
% model.w as parameter. Param is ignored (included for standard
% interface). The loss used is normalized Hamming loss.
% 
% If yi is not given, then standard prediction is done (i.e. MAP decoding
% without the loss term).

if issparse(model.w)
    model.w = full(model.w);
end

% useMex = 1 uses a C-routine for further computation, useMex = 0 uses a   
% Matlab routine (considerably slower!)
useMex = 1; 

%wv = model.w;
num_states = xi.num_states;
num_features = xi.num_features;
featureStart = xi.featureStart;
num_featuresTotal = featureStart(end)-1;
data = xi.data;

if ~useMex
    w = reshape(model.w(1:num_featuresTotal*num_states),num_featuresTotal,num_states);
    v_start = model.w(num_featuresTotal*num_states+1:num_featuresTotal*num_states+num_states);
    v_end = model.w(num_featuresTotal*num_states+num_states+1:num_featuresTotal*num_states+2*num_states);
end
v = reshape(model.w(num_featuresTotal*num_states+2*num_states+1:end),num_states,num_states);

num_Nodes = size(xi.data,1);
    
% Make Potentials
if useMex
    logNodePot = crfChain_makeLogNodePotentialsC(data,model.w,int32(featureStart),int32(num_states));
else
    logNodePot = crfChain_makeLogNodePotentials(xi,w,v_start,v_end,v);
end

% Add loss-augmentation to the score (normalized Hamming distance used for loss)
if nargin > 3
    % Make Loss-Augmented Potentials
    logNodePot = logNodePot + 1/num_Nodes;
    for n = 1:num_Nodes
        logNodePot(n,yi(n)) = logNodePot(n,yi(n)) - 1/num_Nodes;
    end
end
    
% Solve inference problem
yMAP = crfChain_logDecode(logNodePot,v);
label = yMAP;

end