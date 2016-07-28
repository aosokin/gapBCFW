function [components,apps] = modelcomponents(model)
% Function taken from [1]
% [1] Articulated Pose Estimation by a Graphical Model with Image Dependent
% Pairwise Relations, Chen and Yuille

components = cell(length(model.components),1);
for c = 1:length(model.components)
  for k = 1:length(model.components{c})
    p = model.components{c}(k); % has nbh_IDs
    nbh_N = numel(p.nbh_IDs);
    p.Im = cell(nbh_N,1);
    % store the scale of each part relative to the component root
    par = p.parent;
    assert(par < k);
    p.b = [model.bias(p.biasid).w];
    p.b = reshape(p.b,[1 size(p.biasid)]);
    p.biasI = [model.bias(p.biasid).i];
    p.biasI = reshape(p.biasI,size(p.biasid));
    
    x = model.apps(p.appid);
    
    p.sizy = model.tsize(1);
    p.sizx = model.tsize(2);
    p.appI = x.i;
    
    for d = 1:nbh_N
      x = model.pdefs(p.pdefid(d));
      p.pdw(d) = x.w;
      p.pdefI(d) = x.i;
    end
    
    for d = 1:nbh_N
      for m = 1:numel(p.gauid{d})
        x = model.gaus(p.gauid{d}(m));
        p.gauw{d}(m,:)  = x.w;
        p.gauI{d}(m) = x.i;
        mean_x = x.mean(1);
        mean_y = x.mean(2);
        var_x = x.var(1);
        var_y = x.var(2);
        
        p.mean_x{d}(m) = mean_x;
        p.mean_y{d}(m) = mean_y;
        
        p.var_x{d}(m) = var_x;
        p.var_y{d}(m) = var_y;
      end
    end
    components{c}(k) = p;
  end
end
apps = cell(length(model.apps),1);

for i = 1:length(apps)
  apps{i} = model.apps(i).w;
end
