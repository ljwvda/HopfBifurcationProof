function Y = Ybound_G(A,Gext,symmetry,shape,nu,eta,etaPhase)
% Ybound_G computes Y0 bound (residual bound) for existence theorem
%
% Inputs:
%   A - Approximate inverse of Jacobian
%   Gext - Nonlinear function evaluation G
%   symmetry - Symmetry group 
%   shape - Index set shape for bounds
%   nu - Viscosity parameter
%   eta - Weights for x,y
%   etaPhase - Phase condition weights [etaOmega, etaNu]

% Extract phase condition weights
etaOmega = etaPhase(1);
etaNu = etaPhase(2);

% Generate index grids
[nx, ny, nz, nt, ninshape, comp] = countersymtensor(sizeshape_Hopf(shape), shape);

% Handle special case for Hopf bifurcation with nz=2 or nz=1
if strcmp(shape.type, 'ellHopf')
    nz = 2*ones(size(nz));  % Force nz=2 for Hopf case
elseif strcmp(shape.type, 'ellHopf_nz1')
    nz = ones(size(nz));   % Force nz=1 for nz1 Hopf case
end

% Adjust symmetry index for ellHopf/ellHopf_nz1 case to include both vector components
symmetry_idx = symmetry;
if (strcmp(shape.type, 'ellHopf') || strcmp(shape.type, 'ellHopf_nz1')) && symmetry == 1
    symmetry_idx = 11;
elseif (strcmp(shape.type, 'ellHopf') || strcmp(shape.type, 'ellHopf_nz1')) && symmetry == 25
    symmetry_idx = 26;
end

[symvar, symindex, ~, grouporder, multiplicity] = ...
    symmetryindicestensor_ext(nx, ny, nz, nt, comp, sizeshape_Hopf(shape), ninshape, symmetry_idx);

% Identify relevant nz modes (nz=2 for ellHopf, nz=1 for ellHopf_nz1)
if strcmp(shape.type, 'ellHopf_nz1')
    nz_select = (nz(symvar) == 1);
else
    nz_select = (nz(symvar) == 2);
end

% Compute wavenumber squared tensor
[tilden2, ~] = tildentensor(nx, ny, nz, nu);

% Partition indices into finite and tail parts
nfinite = shapetensor(nx, ny, nz, nt, shape); 
finitesymvar = (symvar & nfinite);
symindexjac = symindex(finitesymvar);
tailsymvar = (symvar & ~nfinite);
symindextail = symindex(tailsymvar);  % Used later in function

% Compute weight norm using orbit-stabilizer formula
orbits = grouporder./multiplicity;
weights = eta.^(abs(nx) + abs(ny) + abs(nz) + abs(nt)).*orbits;
weights = weights(symvar);
weights = weights(symindexjac);
finiteweightseta = weights(nz_select);  % Weights for relevant nz modes

% Eigenvalue reciprocals (time-independent case: nt=0)
lambda = reshape(1./abs(nu*tilden2(:)), size(tilden2));


%% Compute weighted residual bound

% Split nonlinear function into real and imaginary parts
Nr = (length(Gext) - 2)/2;
Gextr = Gext(1:Nr);              % Real part (excluding phase equations)
Gexti = Gext(Nr+1:end-2);        % Imaginary part (excluding phase equations)

% Apply approximate inverse to full system
AG = A*Gext;
NAG = (length(AG) - 2)/2;
AGr1 = AG(1:NAG); % Real part of A*G

% Compute Y0 bound: weighted norm of A*G
% Include real part with spatial weights and frequency phase condition
Y1r = sum(abs([AGr1; AG(end-1)]).*[finiteweightseta; etaOmega]);

AGi1 = AG(NAG+1:end-2); % Imaginary part of A*G

% Include imaginary part with spatial weights and viscosity phase condition
Y1i = sum(abs([AGi1; AG(end)]).*[finiteweightseta; etaNu]);

% Tail part contributions (modes beyond finite truncation)
% Real part tail residual
tailresiduer = Gextr(symindextail).*lambda(symindextail);
Y2r = sum(abs(tailresiduer).*weights(symindextail));

% Imaginary part tail residual  
tailresiduei = Gexti(symindextail).*lambda(symindextail);
Y2i = sum(abs(tailresiduei).*weights(symindextail));

% Total Y0 bound: sum of finite and tail contributions
% Note: tail parts are typically zero since G is computed only on finite modes
Y = Y1r + Y1i + Y2r + Y2i;

end