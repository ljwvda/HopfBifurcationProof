function Z2=Z2bound_G(A,rv1,iv1,symmetry,shape,nu,eta,etaPhase,Ndagger,r)
% Z2bound_G  Compute Z2 bound  for Newton-Kantorovich theorem  
%
% Inputs:
%   A         - Approximate inverse of Jacobian  
%   rv1, iv1  - Real and imaginary parts of eigenvector
%   symmetry  - Symmetry group 
%   shape     - Index set shape for bounds
%   nu        - Viscosity parameter
%   eta       - Weights for x,y
%   etaPhase  - Phase condition weights [etaOmega, etaNu]
%   Ndagger   - Finite part truncation parameter
%   r         - Radius for bound computation

% Extract phase condition weights
etaOmega = etaPhase(1);
etaNu = etaPhase(2);

% Generate index grids and symmetry information
M = sizeshape_Hopf(shape);
[nx, ny, nz, nt, ninshape, comp] = countersymtensor(M, shape);

% Handle special case for Hopf bifurcation with nz=2
if strcmp(shape.type, 'ellHopf')
    nz = 2*ones(size(nz));  % Force nz=2 for Hopf case
end

% Adjust symmetry index for ellHopf case
symmetry_idx = symmetry;
if strcmp(shape.type, 'ellHopf') && symmetry == 1
    symmetry_idx = 11;  % Include components 1 and 2 for nz=2
elseif strcmp(shape.type, 'ellHopf') && symmetry == 25
    symmetry_idx = 26;  % Include components 1 and 2 for nz=2
end

[symvar, symindex, ~, grouporder, multiplicity] = ...
    symmetryindicestensor_ext(nx, ny, nz, nt, comp, M, ninshape, symmetry_idx);

% Identify finite part indices
nfinite = shapetensor(nx, ny, nz, nt, shape);
finitesymvar = (symvar & nfinite);
symindexfinite = symindex(finitesymvar);

% Compute weights using orbit-stabilizer formula
orbits = grouporder./multiplicity;
weights = eta.^(abs(nx) + abs(ny) + abs(nz) + abs(nt)).*orbits;
weights = weights(symvar);
finiteweightsetetaPhase = weights(symindexfinite);
nz2 = (nz(symvar)==2); % Select only nz == 2 components
% finiteweightsetetaOmega = finiteweightsetetaOmega(nz2);
finiteweightsetetaPhase = [finiteweightsetetaPhase; etaOmega; etaNu]; 

% operator norm in the rescaled symmetrized norm
% excluding the phase equation
% but including the frequency variable
v=finiteweightsetetaPhase'; 
vv=v(1:end-2);
maxn=max(max(max(max(1,abs(nx)),abs(ny)),abs(nz)),abs(nt));
vvmaxn=vv./maxn(finitesymvar)';

%% First summation

% Subtract the A matrices needed
N_A = (length(A)-2)/2;
A11=A(1:N_A,1:N_A);
A21=A(N_A+1:end-2,1:N_A);
A31=A(end-1,1:N_A);
A41=A(end,1:N_A);

[tilden2,~] = tildentensor(nx,ny,nz,nu);
tilden2mat = tilden2(finitesymvar);
tilden2mat = diag(tilden2mat);

% sqrt(2) for intervals and no intervals
if exist('intval','file') && isintval(nu)
    roottwo=sqrt(intval(2));
else
    roottwo=sqrt(2);
end

normA11maxn=max(v(1:end-2)*abs(A11)./vvmaxn);
normA21maxn=max(v(1:end-2)*abs(A21)./vvmaxn);
normA31=sum(abs(A31)./vv)*etaOmega; 
normA41=sum(abs(A41)./vv)*etaNu; 

Z2sum1_1 = max(v(1:end-2)*abs(A11*tilden2mat)./vvmaxn) + ...
    max(v(1:end-2)*abs(A21*tilden2mat)./vvmaxn) + ...
    sum(abs(A31*tilden2mat)./vv)*etaOmega + ...
    sum(abs(A41*tilden2mat)./vv)*etaNu;
Z2sum1_1=2*Z2sum1_1*(1/etaNu);

Z2sum1_2 = 4*(2/nu)*(1/etaNu); % summation over four terms
normx = sum(abs(rv1(nz2)).*vv');
normy = sum(abs(iv1(nz2)).*vv');
Z2sum1_3aux = (normx+normy+2*r)*(1/nu^3+1/(4*nu^4))*4*(3+roottwo)*eta^2;
% The diagonal is only present in A11, so only there we need the max 
Z2sum1_3 = (4*2*(3+roottwo)*eta^2/nu^2*(1/etaNu)+2*(1/etaOmega)+Z2sum1_3aux*(1/etaNu^2))*(max([normA11maxn,1/(sqrt(nu*Ndagger))])+normA21maxn+normA31+normA41);  

%% Second summation

% Subtract the A matrices needed
A12=A(1:N_A,N_A+1:end-2);
A22=A(N_A+1:end-2,N_A+1:end-2);
A32=A(end-1,N_A+1:end-2);
A42=A(end,N_A+1:end-2);

normA12maxn=max(v(1:end-2)*abs(A12)./vvmaxn);
normA22maxn=max(v(1:end-2)*abs(A22)./vvmaxn);
normA32=sum(abs(A32)./vv)*etaOmega; % No need to use vvmaxn here (as we have vector)
normA42=sum(abs(A42)./vv)*etaNu; % No need to use vvmaxn here


Z2sum2_1 = max(v(1:end-2)*abs(A12*tilden2mat)./vvmaxn) + ...
    max(v(1:end-2)*abs(A22*tilden2mat)./vvmaxn) + ...
    sum(abs(A32*tilden2mat)./vv)*etaOmega + ...
    sum(abs(A42*tilden2mat)./vv)*etaNu;
Z2sum2_1=2*Z2sum2_1*(1/etaNu);

Z2sum2_2 = 2/nu*(1/etaNu); % Same as Z2sum1_2

Z2sum2_3aux = (normx+normy+2*r)*(1/nu^3+1/(4*nu^4))*4*(3+roottwo)*eta^2; % Same as Z2sum1_3
% The diagonal is only present in A22, so only there we need the max 
Z2sum2_3 = (4*2*(3+roottwo)*eta^2/nu^2*(1/etaNu)+2*(1/etaOmega)+Z2sum2_3aux*(1/etaNu^2))*(normA12maxn+max([normA22maxn,1/(sqrt(nu*Ndagger))])+normA32+normA42);  

%% Combining
Z2=Z2sum1_1+Z2sum1_2+Z2sum1_3 + Z2sum2_1+Z2sum2_2+Z2sum2_3;
end



