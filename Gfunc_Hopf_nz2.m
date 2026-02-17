function G = Gfunc_Hopf_nz2(x,y,omega,nu,wTr,wTi,N,symmetry,shape)
% Modified version of Gfunc_Hopf.m for only nz=2 components

%% Equilibrium solution
forcing = classicforcingtensor;
if exist('intval','file') && isintval(nu)
    forcing = intval(forcing);
end
w0 = 1/(2*nu)*forcing; 

solshape.type='rec';
solshape.Nrec=N;

w1=fulltosymmetrytensor_ext(setsizetensor(w0,N),solshape,symmetry);
x0=[w1(:);0];
ref=x0(1:end-1); % reference solution
ref=conjugatesymmetric_ext(ref,solshape,symmetry);

%% Jacobian
% Exclude phase condition and use optimized version for nz==2 only
[~,J]=FDFsymflextensor_ext_Hopf_trunc_nz2(x0,nu,forcing,ref,solshape,symmetry,shape);

%% To take out correct part for nz=2
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(sizeshape_Hopf(shape),shape);
if strcmp(shape.type, 'ellHopf') % Specific for Hopf bifurcation problem with nz=2
      nz = 2*ones(size(nz));
end

if (strcmp(shape.type, 'ellHopf') & symmetry == 1)
    symmetry = 11;  % This ensures we include components 1 and 2 for nz=2
elseif (strcmp(shape.type, 'ellHopf') & symmetry == 25)
    symmetry = 26;  % This ensures we include components 1 and 2 for nz=2
end

[symvar,~] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,sizeshape_Hopf(shape),ninshape,symmetry);

nz2 = (nz(symvar)==2);
% J is already computed only for nz==2 elements, so no need to filter J

x = x(nz2);
y = y(nz2);

wTr = wTr(nz2);
wTi = wTi(nz2);

%%
Jr = real(J);
Ji = imag(J);

%% Compute G

% Real part
Gr = Jr*x-Ji*y+y*omega;

% Imaginary part
Gi = Ji*x+Jr*y-omega*x;

%% "Phase condition"
Cp = abs((wTr+1i*wTi)*(wTr+1i*wTi)'); 

% Pr = wTr*x-wTi*y-1; % Real part
Pr = wTr*x-wTi*y-Cp; % Real part
Pi = wTr*y+wTi*x; % Imaginary part

%% Combining everything

G = [Gr; Gi; Pr; Pi];
end