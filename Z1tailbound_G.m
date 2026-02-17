function [Z1tail,Z1tailcoeff]=Z1tailbound_G(wsol,tensors,solshapeOld,symmetry,nu,eta,etaPhase,Omega,w1der,rv1,iv1,shape,Ntilde,Ndagger,setup)
% Z1tailbound_G  Compute Z1 tail bound for Newton-Kantorovich theorem
%
% Inputs:
%   wsol        - Solution vector in symmetry-reduced form
%   tensors     - Structure containing derivative tensors (Dw, Mw, DMw)
%   solshapeOld - Original solution shape
%   symmetry    - Symmetry group 
%   nu          - Viscosity parameter
%   eta         - Weights for x,y
%   etaPhase    - Weights for phase/amplitude conditions [etaOmega, etaNu]
%   Omega       - Frequency parameter
%   w1der       - Derivative of equilibrium w.r.t. nu
%   rv1, iv1    - Real and imaginary parts of eigenvector (tensor form)
%   shape       - Index set shape for finite part
%   Ntilde      - Tail bound parameter
%   Ndagger     - Finite part truncation parameter
%   setup       - Problem dimension ('2D' or 'not2D')
%
% Outputs:
%   Z1tail      - Total Z1 tail bound
%   Z1tailcoeff - Coefficients used in tail estimate

% Extract derivative tensors
Dw = tensors.Dw;   % Spatial derivatives of solution
Mw = tensors.Mw;   % Nonlinear terms
DMw = tensors.DMw; % Derivatives of nonlinear terms

% Get original dimensions and create modified shape 
N = sizeshape(solshapeOld);
solshapeNew = solshapeOld;
solshapeNew.Nrec = [N(1), N(2), 0, 0];  % Set Nz=0, Nt=0 for equilibrium

% Generate index grids 
[nx, ny, nz, nt] = countersymtensor(solshapeNew.Nrec, solshapeNew);

% Compute exponential decay weights
weights = eta.^(abs(nx) + abs(ny) + abs(nz) + abs(nt));

%% Compute weighted norms for tail bound estimates

% Norms of nonlinear terms Mw (convection)
Mwxi = abs(Mw).*weights;
MWxi1 = Mwxi(:,:,:,:,1);
MWxi2 = Mwxi(:,:,:,:,2);
MWxi3 = Mwxi(:,:,:,:,3);
pMw(1) = sum(MWxi1(:));
pMw(2) = sum(MWxi2(:));
pMw(3) = sum(MWxi3(:));

% Norms of derivatives of nonlinear terms DMw
D1Mw = DMw(:,:,:,:,:,1);
D2Mw = DMw(:,:,:,:,:,2);
D3Mw = DMw(:,:,:,:,:,3);
mDMw(1) = sum(abs(D1Mw(:)).*weights(:));
mDMw(2) = sum(abs(D2Mw(:)).*weights(:));
mDMw(3) = sum(abs(D3Mw(:)).*weights(:));

% Norms of spatial derivatives of solution Dw
D1w = Dw(:,:,:,:,:,1);
D2w = Dw(:,:,:,:,:,2);
D3w = Dw(:,:,:,:,:,3);
mDw(1) = sum(abs(D1w(:)).*weights(:));
mDw(2) = sum(abs(D2w(:)).*weights(:));
mDw(3) = sum(abs(D3w(:)).*weights(:));

% Norms of solution components
w = symmetrytofulltensor_ext(wsol, solshapeNew, symmetry);
wxi = abs(w).*weights;
wxi1 = wxi(:,:,:,:,1);
wxi2 = wxi(:,:,:,:,2);
wxi3 = wxi(:,:,:,:,3);
pw(1) = sum(wxi1(:));
pw(2) = sum(wxi2(:));
pw(3) = sum(wxi3(:));

%% Compute tail bound components

% Square root of tail parameter
if exist('intval', 'file') && isintval(nu)
    rootN = sqrt(intval(Ntilde));
else
    rootN = sqrt(Ntilde);
end

% Apply dimension-dependent rescaling for optimal bounds
if exist('setup', 'var') && strcmp(setup, '2D')
    pMw = pMw/sqrt(nu/2);  % 2D scaling factor
else
    pMw = pMw/sqrt(nu/3);  % 3D scaling factor
end

% Compute tail bound components for each vector component
Z1tailm = max(pMw)/rootN + (3*sum(pw)/2 - pw/2 + mDMw + sum(mDw) - mDw)/Ntilde;

% Z1tail1: Maximum over all vector components 
Z1tail1 = 4*max(Z1tailm);
Z1tailcoeff=4*[max(pMw);max(3*sum(pw)/2-pw/2+mDMw+sum(mDw)-mDw)];

Z1tail2 = 2*Omega/Ntilde;

%% Z1tail3

shapep1.type='rec';
Np1size = sizeshape_Hopf(solshapeOld);
shapep1.Nrec = [Np1size(1)+1 Np1size(2)+1 2 0];

% We take fourth output as we only want the quadratic part of the Jacobian,
% evaluated at derivative of equilibrium, which is w1der
[~,~,~,Dpsiderp1]=FDFsymflextensor_ext_Hopf_trunc_nz2([w1der;0],nu,classicforcingtensor,w1der,solshapeOld,symmetry,shapep1);

symmetry_idx = symmetry;
if (strcmp(shape.type, 'ellHopf') & symmetry == 1)
    symmetry_idx = 11;  % This ensures we include components 1 and 2 for nz=2
elseif (strcmp(shape.type, 'ellHopf') & symmetry == 25)
    symmetry_idx = 26;  % This ensures we include components 1 and 2 for nz=2
end

% New size with one more mode, used for Nz=2 only
shapenz2.type='rec';
shapenz2.Nrec = [Np1size(1)+1 Np1size(2)+1 0 0];

[nx_nz2,ny_nz2,nz_nz2,nt_nz2,ninshape_nz2,comp_nz2] = countersymtensor(sizeshape_Hopf(shapenz2),shapenz2);
if strcmp(shape.type, 'ellHopf') % Specific for Hopf bifurcation problem with nz=2
      nz_nz2 = 2*ones(size(nz_nz2));
end
[symvar_nz2,symindex_nz2,~,grouporder_nz2,multiplicity_nz2] = symmetryindicestensor_ext(nx_nz2,ny_nz2,nz_nz2,nt_nz2,comp_nz2,sizeshape_Hopf(shapenz2),ninshape_nz2,symmetry_idx);

tilden2_nz2 = tildentensor(nx_nz2,ny_nz2,nz_nz2,nu);

rv1tp1 = setsizetensor(rv1,[Np1size(1)+1 Np1size(2)+1 2 0]);
rv1tp1_nz2 = rv1tp1(:,:,5,:,:);
rv1p1 = fulltosymmetrytensor_ext(rv1tp1_nz2,shapenz2,symmetry_idx); 

iv1tp1 = setsizetensor(iv1,[Np1size(1)+1 Np1size(2)+1 2 0]);
iv1tp1_nz2 = iv1tp1(:,:,5,:,:);
iv1p1 = fulltosymmetrytensor_ext(iv1tp1_nz2,shapenz2,symmetry_idx);

% Compute tail contributions from derivative terms
Z1tail3_1 = real(Dpsiderp1)*rv1p1;
Z1tail3_2 = imag(Dpsiderp1)*iv1p1;  % Zero since Dpsiderp1 is real
Z1tail4_1 = imag(Dpsiderp1)*rv1p1;  % Zero since Dpsiderp1 is real  
Z1tail4_2 = real(Dpsiderp1)*iv1p1;

% Sum contributions for tail bounds
Z1tail3sum = abs(Z1tail3_1 - Z1tail3_2);
Z1tail4sum = abs(Z1tail4_1 + Z1tail4_2);

%% Compute weighted tail bounds for nz=2 modes

% Partition indices into finite and tail parts for nz=2 case
nfinite_nz2 = shapetensor(nx_nz2, ny_nz2, nz_nz2, nt_nz2, shape);
tailsymvar_nz2 = (symvar_nz2 & ~nfinite_nz2);
symindextail_nz2 = symindex_nz2(tailsymvar_nz2);

% Compute weights for nz=2 modes using orbit-stabilizer formula
orbits_nz2 = grouporder_nz2./multiplicity_nz2;
weights_nz2 = eta.^(abs(nx_nz2) + abs(ny_nz2) + abs(nz_nz2) + abs(nt_nz2)).*orbits_nz2;
weights_nz2 = weights_nz2(symvar_nz2);

% Eigenvalue reciprocals for nz=2 modes
lambda_nz2 = reshape(1./abs(nu*tilden2_nz2(:)), size(tilden2_nz2));
lambda_nz2 = lambda_nz2(symvar_nz2);

% Compute final tail bound components
Z1tail3 = sum(abs(lambda_nz2(symindextail_nz2).*Z1tail3sum(symindextail_nz2)).*weights_nz2(symindextail_nz2)); 
Z1tail4 = sum(abs(lambda_nz2(symindextail_nz2).*Z1tail4sum(symindextail_nz2)).*weights_nz2(symindextail_nz2)); 

% Report individual tail bound components
disp(['Z1tail1: ', num2str(altsup(Z1tail1))]);
disp(['Z1tail2: ', num2str(altsup(Z1tail2))]);
disp(['Z1tail3: ', num2str(altsup(Z1tail3))]);
disp(['Z1tail4: ', num2str(altsup(Z1tail4))]);

% Total Z1 tail bound: sum of all components
Z1tail = Z1tail1 + Z1tail2 + (1/etaPhase(2))*Z1tail3 + (1/etaPhase(2))*Z1tail4;

end