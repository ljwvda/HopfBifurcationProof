function G = Gfunc_Hopf(x,y,omega,nu,wTr,wTi,N,symmetry)
% Inputs:
%   x        - Real part of eigenvector
%   y        - Imaginary part of eigenvector
%   omega    - Frequency parameter
%   nu       - Viscosity parameter
%   wTr, wTi - Real and imaginary parts of approximation eigenvector
%   N        - Dimensions
%   symmetry - Symmetry group 

%% Equilibrium solution
w0 = 1/(2*nu)*classicforcingtensor; 

solshape.type = 'rec';
solshape.Nrec = N;

w1 = fulltosymmetrytensor_ext(setsizetensor(w0, N), solshape, symmetry);
x0 = [w1(:); 0];
ref = x0(1:end-1);  % Reference solution (equilibrium)
ref = conjugatesymmetric_ext(ref, solshape, symmetry);

%% Jacobian computation
[~, J] = FDFsymflextensor_ext_Hopf(x0, nu, classicforcingtensor, ref, solshape, symmetry);

J = J(1:end-1, 1:end-1);  % Exclude phase condition row/column

% Split Jacobian into real and imaginary parts
Jr = real(J);
Ji = imag(J);

%% Compute nonlinear function G

% Real part: Jr*x - Ji*y + omega*y
Gr = Jr*x - Ji*y + y*omega;

% Imaginary part: Ji*x + Jr*y - omega*x  
Gi = Ji*x + Jr*y - omega*x;

%% Phase condition equations
% Scaling factor for phase condition
Cp = abs((wTr + 1i*wTi)*(wTr + 1i*wTi)');

% Phase condition components
Pr = wTr*x - wTi*y - Cp;  % Real part
Pi = wTr*y + wTi*x;       % Imaginary part

%% Combine all components
G = [Gr; Gi; Pr; Pi];
end