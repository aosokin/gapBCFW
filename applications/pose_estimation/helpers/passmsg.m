function [score,Ix,Iy,Imc,Imp] = passmsg(child, parent, cbid, pbid)
assert(numel(cbid) == 1 && numel(pbid) == 1);
Ny  = size(parent.score,1);
Nx  = size(parent.score,2);

MC = numel(child.gauid{cbid});
MP = numel(parent.gauid{pbid});

[score0,Ix0,Iy0] = deal(zeros(Ny,Nx,MC,MP));
for mc = 1:MC
  for mp = 1:MP
      tmp = double(child.score + (child.pdw(cbid)*child.defMap{cbid}(:,:,mc)));
    [score0(:,:,mc,mp), Ix0(:,:,mc,mp), Iy0(:,:,mc,mp)] = ...
      distance_transform( tmp, ...
      child.gauw{cbid}(mc,:), parent.gauw{pbid}(mp,:), ...
      [child.mean_x{cbid}(mc), child.mean_y{cbid}(mc)], ...
      [child.var_x{cbid}(mc), child.var_y{cbid}(mc)], ...
      [parent.mean_x{pbid}(mp), parent.mean_y{pbid}(mp)], ...
      [parent.var_x{pbid}(mp), parent.var_y{pbid}(mp)], ...
      int32(Nx), int32(Ny) );
    
    score0(:,:,mc,mp) = score0(:,:,mc,mp) + parent.pdw(pbid)*parent.defMap{pbid}(:,:,mp);
  end
end
score = reshape(score0, size(score0, 1),size(score0, 2), MC*MP);
[score, Imcp] = max(score, [], 3);
[Imc, Imp] = ind2sub([MC,MP], Imcp);
[Ix, Iy] = deal(zeros(Ny,Nx));

indI = repmat(1:Ny,1,Nx);
indJ = reshape(repmat(1:Nx,Ny,1),1,Nx*Ny);
Ix(:) = Ix0(sub2ind(size(Ix0),indI',indJ', Imc(:), Imp(:)));
Iy(:) = Iy0(sub2ind(size(Iy0),indI',indJ', Imc(:), Imp(:)));





