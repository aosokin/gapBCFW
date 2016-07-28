function [y] = crfChain_logDecode(logNodePot,logEdgePot)

[num_nodes,num_states] = size(logNodePot);

% Forward Pass
alpha = zeros(num_nodes,num_states);
alpha(1,:) = logNodePot(1,:);
for n = 2:num_nodes 
	tmp = repmat(alpha(n-1,:)',1,num_states) + logEdgePot;
	alpha(n,:) = logNodePot(n,:) + max(tmp);
	[~, mxState(n,:)] = max(tmp);
end

% Backward Pass
y = zeros(num_nodes,1);
[~, y(num_nodes)] = max(alpha(num_nodes,:));
for n = num_nodes-1:-1:1
	y(n) = mxState(n+1,y(n+1));
end

end