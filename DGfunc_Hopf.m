function DG = DGfunc_Hopf(x,y,omega,nu,wTr,wTi,N,symmetry)
% Inputs:
%   x, y     - Real and imaginary parts of eigenvector
%   omega    - Frequency parameter  
%   nu       - Viscosity parameter
%   wTr, wTi - Real and imaginary parts of approximation of eigenvector
%   N        - Dimensions
%   symmetry - Symmetry group 

%% Setup equilibrium solution and shapes
w0 = 1/(2*nu)*classicforcingtensor; 

solshape.type = 'rec';
solshape.Nrec = N;

w1 = fulltosymmetrytensor_ext(setsizetensor(w0, N), solshape, symmetry);
ref = w1(:);  % Reference solution (equilibrium)
ref = conjugatesymmetric_ext(ref, solshape, symmetry);

% Determine largest index set needed
M2N = max(sizeshape_Hopf(solshape), 2*N);

shapelarge.type = 'rec';
shapelarge.Nrec = M2N;
[nx, ny, nz, nt, ninshape, comp] = countersymtensor(M2N, shapelarge);
[symvar, ~] = symmetryindicestensor_ext(nx, ny, nz, nt, comp, M2N, ninshape, symmetry);
[tilden2, ~] = tildentensor(nx, ny, nz, nu);

% Indices for Jacobian computation
njac = shapetensor(nx, ny, nz, nt, solshape);
jacsymvar = (symvar & njac);


%% Compute base Jacobian
[~, J] = FDFsymflextensor_ext_Hopf([w1(:); 0], nu, classicforcingtensor, ref, solshape, symmetry); 
J = J(1:end-1, 1:end-1);  % Exclude phase condition

% Split Jacobian into real and imaginary parts
Jr = real(J); 
Ji = imag(J);

%% Gr (G1)
% Gr = Jr*x-Ji*y+omega*y

% Derivative wrt x 
DG1_x = Jr;
% Derivative wrt y
DG1_y = -Ji+omega*eye(size(J));
% Derivative wrt omega
DG1_om = y;
% Derivative wrt nu 
dudnu = -1/(2*nu^2)*classicforcingtensor;
dudnu = fulltosymmetrytensor_ext(setsizetensor(dudnu,N),solshape,symmetry);

%Set Omega=0 to exclude the time derivative term w.r.t. omega
[~,~,~,Jquad]=FDFsymflextensor_ext([dudnu(:);0],nu,classicforcingtensor,ref,solshape,symmetry);
Jquad = Jquad(1:end-1,1:end-1); 

DJnu_1 = diag(tilden2(jacsymvar)); % First part derivative J wrt nu

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