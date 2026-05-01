function [Z1tail,Z1tailcoeff]=Z1tailbound_G(wsol,tensors,solshapeOld,symmetry,nu,eta,etaPhase,Omega,w1der,rv1,iv1,shape,Ntilde,setup)
% Z1tailbound_G computes Z1 tail bound for existence theorem
%
% Inputs:
%   wsol - Solution vector in symmetry-reduced form
%   tensors - Structure containing derivative tensors (Dw, Mw, DMw)
%   solshapeOld - Original solution shape
%   symmetry - Symmetry group 
%   nu - Viscosity parameter
%   eta - Weights for x,y
%   etaPhase - Weights for phase/amplitude conditions [etaOmega, etaNu]
%   Omega - Frequency parameter
%   w1der - Derivative of equilibrium w.r.t. nu
%   rv1, iv1 - Real and imaginary parts of eigenvector (tensor form)
%   shape - Index set shape for finite part
%   Ntilde - Tail bound parameter
%   setup - Problem dimension ('2D' or 'not2D')
%
% Outputs:
%   Z1tail - Total Z1 tail bound
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

% Determine which nz value to use based on shape type
if strcmp(shape.type, 'ellHopf_nz1')
    hopf_nz = 1;
else
    hopf_nz = 2;  % default for ellHopf
end
shapep1.Nrec = [Np1size(1)+1 Np1size(2)+1 hopf_nz 0];

% We take fourth output as we only want the quadratic part of the Jacobian,
% evaluated at derivative of equilibrium, which is w1der
if strcmp(shape.type, 'ellHopf_nz1')
    [~,~,~,Dpsiderp1]=FDFsymflextensor_ext_Hopf_trunc_nz1([w1der;0],nu,classicforcingtensor,w1der,solshapeOld,symmetry,shapep1);
else
    [~,~,~,Dpsiderp1]=FDFsymflextensor_ext_Hopf_trunc_nz2([w1der;0],nu,classicforcingtensor,w1der,solshapeOld,symmetry,shapep1);
end

symmetry_idx = symmetry;
if (strcmp(shape.type, 'ellHopf') || strcmp(shape.type, 'ellHopf_nz1')) && symmetry == 1
    symmetry_idx = 11;
elseif (strcmp(shape.type, 'ellHopf') || strcmp(shape.type, 'ellHopf_nz1')) && symmetry == 25
    symmetry_idx = 26;
end

% New size with one more mode, for the relevant nz only
shapenz.type='rec';
shapenz.Nrec = [Np1size(1)+1 Np1size(2)+1 0 0];

[nx_nzk,ny_nzk,nz_nzk,nt_nzk,ninshape_nzk,comp_nzk] = countersymtensor(sizeshape_Hopf(shapenz),shapenz);
if strcmp(shape.type, 'ellHopf')       % nz=2 case
    nz_nzk = 2*ones(size(nz_nzk));
elseif strcmp(shape.type, 'ellHopf_nz1') % nz=1 case
    nz_nzk = ones(size(nz_nzk));
end
[symvar_nzk,symindex_nzk,~,grouporder_nzk,multiplicity_nzk] = symmetryindicestensor_ext(nx_nzk,ny_nzk,nz_nzk,nt_nzk,comp_nzk,sizeshape_Hopf(shapenz),ninshape_nzk,symmetry_idx);

tilden2_nzk = tildentensor(nx_nzk,ny_nzk,nz_nzk,nu);

% Extract the nz=hopf_nz Fourier slice: index position = 2*hopf_nz+1 in the nz dimension
nz_slice_idx = 2*hopf_nz + 1;
rv1tp1 = setsizetensor(rv1,[Np1size(1)+1 Np1size(2)+1 hopf_nz 0]);
rv1tp1_nzk = rv1tp1(:,:,nz_slice_idx,:,:);
rv1p1 = fulltosymmetrytensor_ext(rv1tp1_nzk,shapenz,symmetry_idx); 

iv1tp1 = setsizetensor(iv1,[Np1size(1)+1 Np1size(2)+1 hopf_nz 0]);
iv1tp1_nzk = iv1tp1(:,:,nz_slice_idx,:,:);
iv1p1 = fulltosymmetrytensor_ext(iv1tp1_nzk,shapenz,symmetry_idx);

% Compute tail contributions from derivative terms
Z1tail3_1 = real(Dpsiderp1)*rv1p1;
Z1tail3_2 = imag(Dpsiderp1)*iv1p1;  % Zero since Dpsiderp1 is real
Z1tail4_1 = imag(Dpsiderp1)*rv1p1;  % Zero since Dpsiderp1 is real  
Z1tail4_2 = real(Dpsiderp1)*iv1p1;

% Sum contributions for tail bounds
Z1tail3sum = abs(Z1tail3_1 - Z1tail3_2);
Z1tail4sum = abs(Z1tail4_1 + Z1tail4_2);

%% Compute weighted tail bounds for relevant nz modes

% Partition indices into finite and tail parts
nfinite_nzk = shapetensor(nx_nzk, ny_nzk, nz_nzk, nt_nzk, shape);
tailsymvar_nzk = (symvar_nzk & ~nfinite_nzk);
symindextail_nzk = symindex_nzk(tailsymvar_nzk);

% Compute weights using orbit-stabilizer formula
orbits_nzk = grouporder_nzk./multiplicity_nzk;
weights_nzk = eta.^(abs(nx_nzk) + abs(ny_nzk) + abs(nz_nzk) + abs(nt_nzk)).*orbits_nzk;
weights_nzk = weights_nzk(symvar_nzk);

% Eigenvalue reciprocals
lambda_nzk = reshape(1./abs(nu*tilden2_nzk(:)), size(tilden2_nzk));
lambda_nzk = lambda_nzk(symvar_nzk);

% Compute final tail bound components
Z1tail3 = sum(abs(lambda_nzk(symindextail_nzk).*Z1tail3sum(symindextail_nzk)).*weights_nzk(symindextail_nzk)); 
Z1tail4 = sum(abs(lambda_nzk(symindextail_nzk).*Z1tail4sum(symindextail_nzk)).*weights_nzk(symindextail_nzk)); 

% Report individual tail bound components
disp(['Z1tail1: ', num2str(altsup(Z1tail1))]);
disp(['Z1tail2: ', num2str(altsup(Z1tail2))]);
disp(['Z1tail3: ', num2str(altsup(Z1tail3))]);
disp(['Z1tail4: ', num2str(altsup(Z1tail4))]);

% Total Z1 tail bound: sum of all components
Z1tail = Z1tail1 + Z1tail2 + (1/etaPhase(2))*Z1tail3 + (1/etaPhase(2))*Z1tail4;

end