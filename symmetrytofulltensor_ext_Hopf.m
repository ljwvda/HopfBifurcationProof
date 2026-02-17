function wout = symmetrytofulltensor_ext_Hopf(win,shape,symmetry)
% converts a vector of symmetry variables to a full tensor
% it requires the shape of the symmetry variables as an input
% the output tensor then also has size N=sizeshape(shape)

% Adjusted version of symmetrytofulltensor_ext.m for Hopf bifurcaiton
% problem

N=sizeshape_Hopf(shape);

[nx,ny,nz,nt,ninshape,comp] = countersymtensor(N,shape);
if strcmp(shape.type, 'ellHopf') % Specific for Hopf bifurcation problem with nz=2
      nz = 2*ones(size(nz));
end
symmetry_idx = symmetry;
if (strcmp(shape.type, 'ellHopf') & symmetry == 1)
    symmetry_idx = 11;  % This ensures we include components 1 and 2 for nz=2
elseif (strcmp(shape.type, 'ellHopf') & symmetry == 25)
    symmetry_idx = 26;  % This ensures we include components 1 and 2 for nz=2
end
[symvar,~,~,~,multiplicity] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,N,ninshape,symmetry);

wtmp=altzeros([2*N+1,3],win(1)); 
wtmp(symvar)=win;
wtmp=symmetrizetensor_ext(wtmp,symmetry);
wout=reshape(wtmp(:)./multiplicity(:),size(wtmp)); % reshaping to please intlab

end

