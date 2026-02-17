function wout = fulltosymmetrytensor_ext(win,shape,symmetry)
% converts a full tensor to a vector of symmetry variables

sz=size(win);
N=(sz(1:4)-1)/2;

[nx,ny,nz,nt,ninshape,comp]=countersymtensor(N,shape);
symvar=symmetryindicestensor_ext(nx,ny,nz,nt,comp,N,ninshape,symmetry);

wout=win(symvar);
 
end
 