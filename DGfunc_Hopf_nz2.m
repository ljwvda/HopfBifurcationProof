function DG = DGfunc_Hopf_nz2(x,y,omega,nu,wTr,wTi,N,symmetry,shape)
%% Computes the DG matrix for the Hopf bifurcation system (nz=2)
% Inputs:
%   x, y - eigenvectors (not yet filtered for nz=2)
%   omega - frequency
%   nu - viscosity 
%   wTr, wTi - real and imaginary parts of wT (not yet filtered for nz=2)
%   N - truncation parameters [N1, N2, N3, N4]
%   symmetry - symmetry parameter
%   shape - shape structure for DG matrix
%
% Output:
%   DG - DG matrix

solshape.type='rec';
solshape.Nrec=N;

%% Jacobian
forcing = classicforcingtensor;
if exist('intval','file') && isintval(nu)
    forcing = intval(forcing);
end
w0 = 1/(2*nu)*forcing; 

w1=fulltosymmetrytensor_ext(setsizetensor(w0,N),solshape,symmetry);
ref=w1(:); % reference solution
ref=conjugatesymmetric_ext(ref,solshape,symmetry);

% Phase condition is already removed here
[~,J]=FDFsymflextensor_ext_Hopf_trunc_nz2([w1(:);0],nu,forcing,ref,solshape,symmetry,shape); 


%% To take out correct part for nz=2
% biggest index set needed
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(sizeshape_Hopf(shape),shape);
if strcmp(shape.type, 'ellHopf') % Specific for Hopf bifurcation problem with nz=2
      nz = 2*ones(size(nz));
end
symmetry_idx = symmetry;
if (strcmp(shape.type, 'ellHopf') & symmetry == 1)
    symmetry_idx = 11;  % This ensures we include components 1 and 2 for nz=2
elseif (strcmp(shape.type, 'ellHopf') & symmetry == 25)
    symmetry_idx = 26;  % This ensures we include components 1 and 2 for nz=2
end
[symvar,~] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,sizeshape_Hopf(shape),ninshape,symmetry_idx);

nz2 = (nz(symvar)==2);
% J is already computed only for nz==2 elements by the optimized function, so no need to filter J

x = x(nz2);
y = y(nz2);

wTr = wTr(nz2);
wTi = wTi(nz2);

%%

Jr = real(J); Ji = imag(J);

%% Gr (G1)
% Gr = Jr*x-Ji*y+omega*y

% Derivative wrt x 
DG1_x = Jr;
% Derivative wrt y
DG1_y = -Ji+omega*eye(size(J));
% Derivative wrt omega
DG1_om = y;
% Derivative wrt nu 
dudnu = -1/(2*nu^2)*forcing;
dudnu = fulltosymmetrytensor_ext(setsizetensor(dudnu,N),solshape,symmetry);

%Set Omega=0 to exclude the time derivative term w.r.t. omega
%%
[~,~,~,Jquad]=FDFsymflextensor_ext_Hopf_trunc_nz2([dudnu(:);0],nu,forcing,ref,solshape,symmetry,shape);

[tilden2,~] = tildentensor(nx,ny,nz,nu);

% indices for the jacobian - filter for nz==2
njac = shapetensor(nx,ny,nz,nt,shape);
jacsymvar = (symvar & njac & (nz==2));

DJnu_1 = diag(tilden2(jacsymvar)); % First part derivative J wrt nu

% DJnu_1 is already filtered for nz==2, so no additional filtering needed

DJr_nu = real(DJnu_1)+real(Jquad); % Derivative real(J) wrt nu
DJi_nu = imag(DJnu_1)+imag(Jquad); % Derivative imag(J) wrt nu

DG1_nu = DJr_nu*x-DJi_nu*y;

DG1 = [DG1_x DG1_y DG1_om DG1_nu];

%% Gi (G2)
% Gi = Ji*x+Jr*y-omega*x

% Derivative wrt x
DG2_x = Ji-omega*eye(size(J));
% Derivative wrt y
DG2_y = Jr;
% Derivative wrt omega
DG2_om = -x;
% Derivative wrt nu
DG2_nu = DJi_nu*x+DJr_nu*y;

DG2 = [DG2_x DG2_y DG2_om DG2_nu];

%% Pr (G3)
% Pr = real(wT)*x-imag(wT)*y-1

% Derivative wrt x
DG3_x = wTr;
% Derivative wrt y 
DG3_y = -wTi;
% Derivative wrt omega
DG3_om = 0;
% Derivative wrt nu
DG3_nu = 0;

DG3 = [DG3_x DG3_y DG3_om DG3_nu];

%% Pi (G4)
% Pi = real(wT)*y+imag(wT)*x

% Derivative wrt x
DG4_x = wTi;
% Derivative wrt y 
DG4_y = wTr;
% Derivative wrt omega
DG4_om = 0;
% Derivative wrt nu
DG4_nu = 0;

DG4 = [DG4_x DG4_y DG4_om DG4_nu];

%% Combining

DG = [DG1; DG2; DG3; DG4]; 
end