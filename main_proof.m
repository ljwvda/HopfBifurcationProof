%% Method for proving Hopf bifurcation point in Navier-Stokes equations
tic

use_intervals = false;  % Set true to use interval arithmetic for rigorous bounds

if use_intervals
    disp('Proof is done using intervals')
else
    disp('Proof is done without intervals')
end

% Load initial guess and parameter values
load('eigvec1_thirdHopf_32modes.mat','eigvec1');
v1 = eigvec1;
Omega = 5.013;  % Imaginary part of eigenvalue at Hopf bifurcation
nu = 0.12;      % Viscosity parameter

% Extract real and imaginary parts of eigenvector
rv1 = real(v1);
iv1 = imag(v1);

% Initial tensor dimensions (full size for loading eigenvector)
N_initial = [32, 32, 2, 0];  % Modes: nonzero at nz=2 or nz=-2 only
solshape.type = 'rec';
solshape.Nrec = N_initial;

% Convert eigenvector to tensor format
rv1t = symmetrytofulltensor_ext(rv1, solshape, 1);
iv1t = symmetrytofulltensor_ext(iv1, solshape, 1);

Nbig = 10;  % Truncated modes for faster computation
N = [Nbig, Nbig, 2, 0];
solshape.type = 'rec';
solshape.Nrec = N;

symmetry = 25;  % Including symmetries

% Resize eigenvectors to adjusted dimensions
rv1 = fulltosymmetrytensor_ext(setsizetensor(rv1t, N), solshape, symmetry);
rv1 = rv1(:);

iv1 = fulltosymmetrytensor_ext(setsizetensor(iv1t, N), solshape, symmetry);
iv1 = iv1(:);

% Adjoint eigenvector components (complex conjugate transpose)
wTr = rv1.';   % Real part of adjoint eigenvector
wTi = -iv1.';  % Imaginary part (negative due to complex conjugation)

%% Newton iterations
disp('First run Newton iterations');
tol=1; step=0;
while tol>1e-14 && step<10
    G = Gfunc_Hopf(rv1,iv1,Omega,nu,wTr,wTi,N,symmetry);
    DG = DGfunc_Hopf(rv1,iv1,Omega,nu,wTr,wTi,N,symmetry);
    dx=-DG\G;
    tol=norm(dx);
    disp(['Stepsize is ',num2str(tol)])
    x0 = [rv1;iv1;Omega;nu];
    x0=x0+dx;
    rv1 = x0(1:length(rv1));
    iv1 = x0(length(rv1)+1:2*length(rv1));
    Omega = x0(end-1);
    nu = x0(end);
    step=step+1;
end
disp('End of Newton iterations');

%% Convert to intervals for rigorous computation after Newton convergence
if use_intervals && exist('intval','file')
    rv1 = intval(rv1);
    iv1 = intval(iv1);
    Omega = intval(Omega);
    nu = intval(nu);
end

newv1 = rv1+1i*iv1;

rv1t = symmetrytofulltensor_ext(rv1,solshape,symmetry);
iv1t = symmetrytofulltensor_ext(iv1,solshape,symmetry);

wTr = rv1.';% Complex conjugate
wTi = -iv1.'; % - because we take complex conjugate

% Change N
N=[Nbig,Nbig,0,0]; % nonzero elements are at nz=2 and nz=-2
solshape.type='rec';
solshape.Nrec=N;

%% Preparing for the proof

% Compute equilibrium solution
forcing = classicforcingtensor;
if use_intervals && exist('intval','file')
    forcing = intval(forcing);
end
w0 = 1/(2*nu)*forcing;  % Equilibrium solution of Navier-Stokes
w1 = fulltosymmetrytensor_ext(setsizetensor(w0, N), solshape, symmetry);
w1 = w1(:);

% Weight parameters for norm estimates
eta = 1;         % Spatial decay parameter
etaOmega = 1;    % Weight for frequency parameter
etaNu = 70;      % Weight for viscosity parameter  
etaPhase = [etaOmega; etaNu];

% Set up shape with Z component for tensor computations
solshapeZ.type = 'rec';
NZ = [Nbig, Nbig, 2, 0];
solshapeZ.Nrec = NZ;

[~,~,~,~,~,~,tensors]=FDFsymflextensor_ext_Hopf([w1;Omega],nu,forcing,w1,solshape,symmetry);

setup = '2D'; % As equilibrium is 2D, we can use sqrt(2) instead of sqrt(3)
Z1target = 0.9175;

tildeNestimate=estimatetildeN_G(Z1target,w1,tensors,solshapeZ,symmetry,nu,eta,1,setup); % Set Ntilde=1 for estimate
disp(['The required Ntilde for the tail term is about ',int2str(tildeNestimate)]);

%  Ntilde >= Ndagger (to ensure finite part A does not influence tail
%  estimate)
Ndagger = 20; % Reduced for faster computation
Ntilde = 30; % Reduced for faster computation

%% Create Edaggershape and Etildeshape, specifically for 2D case
Edaggersh='ell2D';
Etildesh='ell2D';

Edaggershape.type=Edaggersh;
Edaggershape.Nell=Ndagger;
Edaggershape.nu=nu;
Edaggershape.omega=Omega;

Etildeshape.type=Etildesh;
Etildeshape.Nell=Ntilde;
Etildeshape.nu=nu;
Etildeshape.omega=Omega;

%% Create shape EdaggershapeHopf
EdaggershHopf='ellHopf';

EdaggershapeHopf.type=EdaggershHopf;
EdaggershapeHopf.Nell=Ndagger;
EdaggershapeHopf.nu=nu;
EdaggershapeHopf.omega=Omega;

%% Make full 3D ellips
Edaggersh3D='ell';

N=sizeshape_Hopf(solshapeZ);
Edaggershape3D.type=Edaggersh3D;
Edaggershape3D.Nell=Ndagger;
Edaggershape3D.nu=nu;
Edaggershape3D.omega=Omega;

%% Preparation for computing G and DG
w0 = 1/(2*nu)*classicforcingtensor;
w1=fulltosymmetrytensor_ext(setsizetensor(w0,N),solshape,symmetry);
ref=w1(:); % reference solution
ref=conjugatesymmetric_ext(ref,solshape,symmetry);

% Change size eigenvectors
rv1t_for_ellips = setsizetensor(rv1t, sizeshape_Hopf(solshapeZ));
iv1t_for_ellips = setsizetensor(iv1t, sizeshape_Hopf(solshapeZ));

% For the ellips case, use Edaggershape3D
rv1tlarge_ell3D = setsizetensor(rv1t_for_ellips, sizeshape_Hopf(Edaggershape3D));
iv1tlarge_ell3D = setsizetensor(iv1t_for_ellips, sizeshape_Hopf(Edaggershape3D));

% Go to vector format using the 3D ellips shape
rv1large_ell3D = fulltosymmetrytensor_ext(rv1tlarge_ell3D, Edaggershape3D, symmetry);
iv1large_ell3D = fulltosymmetrytensor_ext(iv1tlarge_ell3D, Edaggershape3D, symmetry);

wTrlarge_ell3D = rv1large_ell3D.';% Complex conjugate
wTilarge_ell3D = -iv1large_ell3D.'; % - because we take complex conjugate

% The function will filter to nz=2 components and use the 3D ellips shape for Jacobian
Glarge = Gfunc_Hopf_nz2(rv1large_ell3D,iv1large_ell3D,Omega,nu,wTrlarge_ell3D,wTilarge_ell3D,sizeshape_Hopf(Edaggershape3D),symmetry,Edaggershape3D);
disp(['Norm of Glarge (Edaggershape3D) is ',num2str(altsup(norm(Glarge)))]);

DGext = DGfunc_Hopf_nz2(rv1large_ell3D,iv1large_ell3D,Omega,nu,wTrlarge_ell3D,wTilarge_ell3D,sizeshape_Hopf(Edaggershape3D),symmetry,Edaggershape3D);

%% Calculate A (convert to intervals if requested)
if use_intervals && exist('intval','file')
    Aext=intval(inv(mid(DGext)));
else
    Aext=inv(DGext);
   % disp(['norm of A is ',num2str(norm(Aext,1))]);
end

%% Bounds
%%% Y0 bound %%%
disp('Start computation Y0 bound');
Y0 = Ybound_G(Aext,Glarge,symmetry,EdaggershapeHopf,nu,eta,etaPhase);
disp(['Y0 is ',num2str(altsup(Y0))]);

%%% Z1 bound %%%
disp('Start computation Z1 bound');

%%% Z1tail %%%%
w0der = -1/(2*nu^2)*forcing;  % derivative of equilibrium solution
w1der=fulltosymmetrytensor_ext(setsizetensor(w0der,N),solshape,symmetry);
w1der=w1der(:);

Z1tail=Z1tailbound_G(w1,tensors,solshape,symmetry,nu,eta,etaPhase,Omega,w1der,rv1t,iv1t,EdaggershapeHopf,Ntilde,Ndagger,setup);
disp(['the Z1 tail term is ',num2str(altsup(Z1tail))]);

%%% Z2 %%%
if use_intervals
    rtest = intval('0.001');
else
    rtest = 0.001; 
end
Z2test=Z2bound_G(Aext,rv1large_ell3D,iv1large_ell3D,symmetry,EdaggershapeHopf,nu,eta,etaPhase,Ndagger,rtest);
X = ['The value of the test Z2-bound is ',num2str(altsup(Z2test)), ' for rtest = ', num2str(altsup(rtest)),'.'];
disp(X)

%%% Z1finite %%%
% Create Q2Nshape and QNshape for the new DG computation
QNshape.type='otherplusrec'; % Number of columns
QNshape.other=Etildeshape;
QNshape.Nrec=N;

Q2Nshape.type='otherplusrec'; % Number of rows
Q2Nshape.other=Etildeshape;
Q2Nshape.Nrec=2*N;

% Create vectors compatible with Q2Nshape for DG computation
% We need to create rv1, iv1, wTr, wTi that are compatible with Q2Nshape
rv1tlarge_Q2N = setsizetensor(rv1t_for_ellips, sizeshape_Hopf(Q2Nshape));
iv1tlarge_Q2N = setsizetensor(iv1t_for_ellips, sizeshape_Hopf(Q2Nshape));

% Go to vector format using the Q2Nshape
rv1large_Q2N = fulltosymmetrytensor_ext(rv1tlarge_Q2N, Q2Nshape, symmetry);
iv1large_Q2N = fulltosymmetrytensor_ext(iv1tlarge_Q2N, Q2Nshape, symmetry);

wTrlarge_Q2N = rv1large_Q2N.';% Complex conjugate
wTilarge_Q2N = -iv1large_Q2N.'; % - because we take complex conjugate

%% Parallel version
Z1finite_parallel = Z1finitebound_G_parallel(Aext,rv1large_Q2N,iv1large_Q2N,Omega,nu,wTrlarge_Q2N,wTilarge_Q2N,symmetry,Edaggershape3D,N,Etildeshape,eta,etaPhase);
disp(['Z1 finite is (parallel) ',num2str(altsup(Z1finite_parallel))]);
Z1finite = Z1finite_parallel;
% Z1finite = Z1finitebound_G(Aext,rv1large_Q2N,iv1large_Q2N,Omega,nu,wTrlarge_Q2N,wTilarge_Q2N,symmetry,Edaggershape3D,N,Etildeshape,eta,etaPhase);
% disp(['Z1 finite is ',num2str(altsup(Z1finite))]);

Z1=max([Z1finite,Z1tail]);

%% Start the proof

% Convert bounds to floating points
Y0 = altsup(Y0);
Z1 = altsup(Z1);
Z2test = altsup(Z2test);

discr_test=(1-Z1)^2-2*Y0*Z2test;
r0 = altsup((1-Z1-sqrt(discr_test))/Z2test);
% We can compute Z2(r0) once r0 is fixed
Z2=Z2bound_G(Aext,rv1large_ell3D,iv1large_ell3D,symmetry,EdaggershapeHopf,nu,eta,etaPhase,Ndagger,r0);
X = ['The value of the Z2-bound is ',num2str(altsup(Z2)), ' for r0 = ', num2str(altsup(r0)),'.'];
disp(X)

discr=(1-Z1)^2-2*Y0*Z2;
disp(['discriminant is ',num2str(altsup(discr))]);
% Start the actual proof
if 1-Z1>0 && discr>0
   disp('SUCCESS')
   if exist('intval','file') && isintval(nu)
      disp('including full interval arithmetic')
   else
      disp('but no interval arithmetic')
   end
   rmin=altsup((1-Z1-sqrt(discr))/Z2);
   rmax=-altsup(-(1-Z1)/Z2);
   disp(['r_min is ',num2str(rmin)]);
   disp(['r_max is ',num2str(rmax)]);
   success=true;
else
   disp('FAILURE')
end

toc