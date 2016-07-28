function phi = chain_conll_featuremap(param, xi, yi)
% computes the joint feature map phi(x,y). [param is ignored]

%Determines whether to use C routine for the featuremap computation (this
%are kept to see what is happening in the C code)
useMex = 1; 

num_states = xi.num_states;

featureStart = xi.featureStart;
num_featuresTotal = featureStart(end)-1;
data = xi.data;

unit = zeros(num_featuresTotal,num_states);
gr_start = zeros(num_states,1);
gr_end = zeros(num_states,1);
bin = zeros(num_states);

% Note that xi.data is given already by its unary features. In the 
% following, these are transformed into a binary representation and biases
% for the beginning and the end of the sentence as well as sequential 
% features are added.

if useMex
    crfChain_maxMarginUpdateGradient(unit,gr_start,gr_end,bin,data,...
        int32(yi),num_states,featureStart);
else
    num_features = xi.num_features;
    nNodes = size(xi.data,1);
    % Update gradient
    for n = 1:nNodes % unary features
        features = xi.data(n,:);
        for feat = 1:length(num_features)
            if features(feat) ~= 0
                featureParam = featureStart(feat)+features(feat)-1;
                for state = 1:num_states
                    O = (state == yi(n)); 
                    unit(featureParam,state) = unit(featureParam,state) + O;
                end
             end
         end
    end
    for state = 1:num_states
        O = (state == yi(1));
        gr_start(state) = gr_start(state) + O; % beginning of sentence
        O = (state == yi(end));
        gr_end(state) = gr_end(state) + O; % end of sentence
    end
    for n = 1:nNodes-1 % sequential binary features
        for state1 = 1:num_states
            for state2 = 1:num_states
                O = ((state1 == yi(n)) && (state2 == yi(n+1)));
                bin(state1,state2) = bin(state1,state2) + O;
            end
        end
    end
end

phi = [unit(:);gr_start;gr_end;bin(:)];
% transform into space feature
phi = sparse(phi);

end
