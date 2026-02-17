function [nsymvar,NN]=makeindices(N,nx,ny,nz,nt,symvar)
% constructs some Fourier indices and sizes
% used in the efficient computation of derivative terms
% in particular for turning tensor indices into "linear" indices

nsymvar=zeros([length(find(symvar)),4]);
nsymvar(:,1)=nx(symvar);
nsymvar(:,2)=ny(symvar);
nsymvar(:,3)=nz(symvar); 
nsymvar(:,4)=nt(symvar);

NN=zeros([5,1]);
NN(1)=2*N(1)+1;
NN(2)=NN(1)*(2*N(2)+1);
NN(3)=NN(2)*(2*N(3)+1);
NN(4)=NN(3)*(2*N(4)+1);
NN(5)=3*NN(4);

end
