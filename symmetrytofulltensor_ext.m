function wout = symmetrytofulltensor_ext(win,shape,symmetry)
% converts a vector of symmetry variables to a full tensor
% it requires the shape of the symmetry variables as an input
% the output tensor then also has size N=sizeshape(shape)

N=sizeshape(shape);

[nx,ny,nz,nt,ninshape,comp] = countersymtensor(N,shape);
[symvar,~,~,~,multiplicity] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,N,ninshape,symmetry);

wtmp=altzeros([2*N+1,3],win(1)); 
wtmp(symvar)=win;
wtmp=symmetrizetensor_ext(wtmp,symmetry);
wout=reshape(wtmp(:)./multiplicity(:),size(wtmp)); % reshaping to please intlab

end

