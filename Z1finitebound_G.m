function Z1finite = Z1finitebound_G(A,rv1large_Q2N,iv1large_Q2N,Omega,nu,wTrlarge_Q2N,wTilarge_Q2N,symmetry,Edaggershape,N,Etildeshape,eta,etaPhase)
% Z1finitebound_G  Compute Z1 finite bound for Newton-Kantorovich theorem
%
% Inputs:
%   A                    - Approximate inverse of Jacobian
%   rv1large_Q2N, iv1large_Q2N - Extended eigenvectors (real/imaginary)
%   Omega                - Frequency parameter
%   nu                   - Viscosity parameter
%   wTrlarge_Q2N, wTilarge_Q2N - Extended approximation of eigenvectors
%   sizeshape_dagger     - Size of finite part shape
%   symmetry             - Symmetry group
%   Edaggershape         - 3D ellips shape structure based on Ndagger
%   N                    - Original tensor dimensions
%   Etildeshape          - 3D ellips shape structure based on Ntilde
%   eta                  - Weights for x,y
%   etaPhase             - Phase condition weights [etaOmega, etaNu]

% Set up extended shapes for finite bound computation
QNshape.type = 'otherplusrec';    % Column space
QNshape.other = Etildeshape;
QNshape.Nrec = N;

Q2Nshape.type = 'otherplusrec';   % Row space
Q2Nshape.other = Etildeshape;
Q2Nshape.Nrec = 2*N;

% Extract phase condition weights
etaOmega = etaPhase(1);
etaNu = etaPhase(2);

%% Compute weights and index mappings for different spaces

% Row space (Q2Nshape)
[nx_dg_row, ny_dg_row, nz_dg_row, nt_dg_row, ninshape_dg_row, comp_dg_row] = countersymtensor(sizeshape_Hopf(Q2Nshape), Q2Nshape);
[symvar_dg_row,~,~,~,~,~] = symmetryindicestensor_ext(nx_dg_row, ny_dg_row, nz_dg_row, nt_dg_row, comp_dg_row, sizeshape_Hopf(Q2Nshape), ninshape_dg_row, symmetry);
nz2_dg_row = (nz_dg_row(symvar_dg_row) == 2);

% Column space (QNshape)
[nx_dg_col, ny_dg_col, nz_dg_col, nt_dg_col, ninshape_dg_col, comp_dg_col] = countersymtensor(sizeshape_Hopf(QNshape), QNshape);
[symvar_dg_col, ~, ~, grouporder_dg_col, multiplicity_dg_col, ~] = symmetryindicestensor_ext(nx_dg_col, ny_dg_col, nz_dg_col, nt_dg_col, comp_dg_col, sizeshape_Hopf(QNshape), ninshape_dg_col, symmetry);
nz2_dg_col = (nz_dg_col(symvar_dg_col) == 2);

% Compute weights for column variables (QNshape space)
orbits_dg_col = grouporder_dg_col./multiplicity_dg_col;
weights_dg_col = eta.^(abs(nx_dg_col) + abs(ny_dg_col) + abs(nz_dg_col) + abs(nt_dg_col)).*orbits_dg_col;
weights_nz2_col = weights_dg_col(symvar_dg_col);
weights_nz2_col_only = weights_nz2_col(nz2_dg_col);  % Extract nz=2 components only

% Set up column weights 
% Matrix structure: rows from Q2Nshape, columns from QNshape
% Columns: [real_part; imaginary_part] for nz=2 modes
weightsDG_col = [weights_nz2_col_only; weights_nz2_col_only];

%% Create mappings between Edaggershape and QNshape/Q2Nshape

% Get Edaggershape indices for nz=2 components
[nx_edag,ny_edag,nz_edag,nt_edag,ninshape_edag,comp_edag] = countersymtensor(sizeshape_Hopf(Edaggershape),Edaggershape);
[symvar_edag,~,~,grouporder_edag,multiplicity_edag,~] = symmetryindicestensor_ext(nx_edag,ny_edag,nz_edag,nt_edag,comp_edag,sizeshape_Hopf(Edaggershape),ninshape_edag,symmetry);
nz2_edag = (nz_edag(symvar_edag)==2);

% Extract Edaggershape nz=2 mode indices
nx_edag_nz2 = nx_edag(symvar_edag);
ny_edag_nz2 = ny_edag(symvar_edag);
nz_edag_nz2 = nz_edag(symvar_edag);
nt_edag_nz2 = nt_edag(symvar_edag);
comp_edag_nz2 = comp_edag(symvar_edag);

nx_edag_nz2 = nx_edag_nz2(nz2_edag);
ny_edag_nz2 = ny_edag_nz2(nz2_edag);
nz_edag_nz2 = nz_edag_nz2(nz2_edag);
nt_edag_nz2 = nt_edag_nz2(nz2_edag);
comp_edag_nz2 = comp_edag_nz2(nz2_edag);

% Extract Q2Nshape nz=2 mode indices
nx_Q2N_nz2 = nx_dg_row(symvar_dg_row);
ny_Q2N_nz2 = ny_dg_row(symvar_dg_row);
nz_Q2N_nz2 = nz_dg_row(symvar_dg_row);
nt_Q2N_nz2 = nt_dg_row(symvar_dg_row);
comp_Q2N_nz2 = comp_dg_row(symvar_dg_row);

nx_Q2N_nz2 = nx_Q2N_nz2(nz2_dg_row);
ny_Q2N_nz2 = ny_Q2N_nz2(nz2_dg_row);
nz_Q2N_nz2 = nz_Q2N_nz2(nz2_dg_row);
nt_Q2N_nz2 = nt_Q2N_nz2(nz2_dg_row);
comp_Q2N_nz2 = comp_Q2N_nz2(nz2_dg_row);

% Map each Edaggershape mode to its position in Q2Nshape (-1 if not found)
edag_to_Q2N_mapping = [];
for i = 1:length(nx_edag_nz2)
    match_idx = find(nx_Q2N_nz2 == nx_edag_nz2(i) & ...
        ny_Q2N_nz2 == ny_edag_nz2(i) & ...
        nz_Q2N_nz2 == nz_edag_nz2(i) & ...
        nt_Q2N_nz2 == nt_edag_nz2(i) & ...
        comp_Q2N_nz2 == comp_edag_nz2(i));
    if ~isempty(match_idx)
        edag_to_Q2N_mapping = [edag_to_Q2N_mapping; match_idx]; 
    else
        edag_to_Q2N_mapping = [edag_to_Q2N_mapping; -1];
    end
end

%% Map Edaggershape to QNshape (column space)

n_vars_col = sum(nz2_dg_col); % Number of nz=2 column variables in QNshape

% Extract QNshape nz=2 mode indices
nx_QN_col_nz2 = nx_dg_col(symvar_dg_col);
ny_QN_col_nz2 = ny_dg_col(symvar_dg_col);
nz_QN_col_nz2 = nz_dg_col(symvar_dg_col);
nt_QN_col_nz2 = nt_dg_col(symvar_dg_col);
comp_QN_col_nz2 = comp_dg_col(symvar_dg_col);

nx_QN_col_nz2 = nx_QN_col_nz2(nz2_dg_col);
ny_QN_col_nz2 = ny_QN_col_nz2(nz2_dg_col);
nz_QN_col_nz2 = nz_QN_col_nz2(nz2_dg_col);
nt_QN_col_nz2 = nt_QN_col_nz2(nz2_dg_col);
comp_QN_col_nz2 = comp_QN_col_nz2(nz2_dg_col);

% Map each Edaggershape mode to its position in QNshape (-1 if not found)
edag_to_QN_col_mapping = [];
for i = 1:length(nx_edag_nz2)
    match_idx = find(nx_QN_col_nz2 == nx_edag_nz2(i) & ...
        ny_QN_col_nz2 == ny_edag_nz2(i) & ...
        nz_QN_col_nz2 == nz_edag_nz2(i) & ...
        nt_QN_col_nz2 == nt_edag_nz2(i) & ...
        comp_QN_col_nz2 == comp_edag_nz2(i));
    if ~isempty(match_idx)
        edag_to_QN_col_mapping = [edag_to_QN_col_mapping; match_idx]; 
    else
        edag_to_QN_col_mapping = [edag_to_QN_col_mapping; -1];
    end
end

%% Map Edaggershape to Q2Nshape (row space)

n_vars_row = sum(nz2_dg_row);

% Filter to valid modes (exist in both Edaggershape and Q2Nshape)
row_mapping_mask = (edag_to_Q2N_mapping > 0);
valid_edag_to_Q2N_mapping = edag_to_Q2N_mapping(row_mapping_mask);

% Create index vector for extracting Edaggershape part from DG: [x_edag; y_edag; omega; nu]
edag_row_indices = [valid_edag_to_Q2N_mapping; ...
    valid_edag_to_Q2N_mapping + n_vars_row; ...
    2*n_vars_row+1; ...
    2*n_vars_row+2];

% Weights for Edaggershape structure (matching A matrix)
orbits_edag=grouporder_edag./multiplicity_edag;
weights_edag=eta.^(abs(nx_edag)+abs(ny_edag)+abs(nz_edag)+abs(nt_edag)).*orbits_edag;
weights_nz2_edag = weights_edag(symvar_edag);
weights_nz2_edag_only = weights_nz2_edag(nz2_edag);
weightsDG_edag = [weights_nz2_edag_only; weights_nz2_edag_only; etaOmega; etaNu]; % [x; y; omega; nu] for Edaggershape

%% Precompute matrix dimensions and remainder indices

A_size = size(A,1); % Size of A matrix (from Edaggershape)
DG_rows = 2*sum(nz2_dg_row) + 2; % Size of DG rows (from Q2Nshape): 2*Q2N+2
DG_cols = 2*sum(nz2_dg_col); % Size of DG columns (from QNshape): 2*QN

disp(['A matrix size: ', num2str(A_size), 'x', num2str(A_size)]);
disp(['DG matrix size: ', num2str(DG_rows), 'x', num2str(DG_cols)]);

% Initialize maximum values
max_col_sum = 0;
max_remainder_eigvec_cols = 0; 

% Precompute remainder indices and properties
all_indices = (1:DG_rows)';
remainder_indices = setdiff(all_indices, edag_row_indices);

% Initialize
remainder_nx = [];
remainder_ny = [];
remainder_nz = [];
remainder_nt = [];
remainder_comp = [];

% Map remainder_indices back to the Q2Nshape structure to get (nx,ny,nz,nt)
nz2_indices = find(nz2_dg_row); % Indices of nz=2 components in Q2Nshape

for i = 1:length(remainder_indices)
    idx = remainder_indices(i);
    if idx <= n_vars_row
        % Real part of nz=2 variable 
        if idx <= length(nz2_indices)
            mode_idx = nz2_indices(idx); % Get the actual mode index
            remainder_nx = [remainder_nx; nx_dg_row(mode_idx)];
            remainder_ny = [remainder_ny; ny_dg_row(mode_idx)];
            remainder_nz = [remainder_nz; nz_dg_row(mode_idx)];
            remainder_nt = [remainder_nt; nt_dg_row(mode_idx)];
            remainder_comp = [remainder_comp; comp_dg_row(mode_idx)];
        end
    elseif idx <= 2*n_vars_row
        % Imaginary part of nz=2 variable
        real_idx = idx - n_vars_row;
        if real_idx <= length(nz2_indices)
            mode_idx = nz2_indices(real_idx); % Get the actual mode index
            remainder_nx = [remainder_nx; nx_dg_row(mode_idx)];
            remainder_ny = [remainder_ny; ny_dg_row(mode_idx)];
            remainder_nz = [remainder_nz; nz_dg_row(mode_idx)];
            remainder_nt = [remainder_nt; nt_dg_row(mode_idx)];
            remainder_comp = [remainder_comp; comp_dg_row(mode_idx)];
        end
    end
    % Skip omega and nu terms (indices 2*n_vars_row+1 and 2*n_vars_row+2)
end

% Compute remainder lambda and weights
tilden2_remainder = remainder_nx.^2 + remainder_ny.^2 + remainder_nz.^2 + remainder_nt.^2;
lambda_remainder = 1./abs(nu*tilden2_remainder);
lambda_remainder(nu*tilden2_remainder==0)=0;
remainder_weights = eta.^(abs(remainder_nx)+abs(remainder_ny)+abs(remainder_nz)+abs(remainder_nt));

%% Loop over QNshape columns to compute I - A*DG column by column

% Loop over Stilde(=QNshape) columns
for col_idx = 1:DG_cols
    
    if mod(col_idx, 1000) == 0
        disp(['Computing QNshape column ', num2str(col_idx), ' out of ', num2str(DG_cols), ', current max: ', num2str(max_col_sum)]);
    end

    % Compute single column of DG matrix with Q2Nshape structure
    DG_col = DGfunc_Hopf_nz2_column(rv1large_Q2N,iv1large_Q2N,Omega,nu,wTrlarge_Q2N,wTilarge_Q2N,N,symmetry,Q2Nshape,QNshape,col_idx);

    % Extract Edaggershape part from DG_col to match A matrix size
    DG_col_edag_extracted = DG_col(edag_row_indices); % Extract Edaggershape part
    
    % Create the full-size DG_col_edag vector with zeros for missing modes
    % A matrix expects size [2*n_edag + 2, 1] where n_edag = number of Edaggershape nz=2 modes
    n_edag = sum(nz2_edag);  % Total number of Edaggershape nz=2 modes
    DG_col_edag = altzeros([2*n_edag + 2, 1], nu);
    
    % Fill in the valid parts (modes that exist in both Edaggershape and Q2Nshape)
    % This is necessary for the case Ndagger >> Ntilde
    n_valid_modes = sum(row_mapping_mask);
    valid_edag_indices = find(row_mapping_mask); % Edaggershape mode indices that exist in Q2Nshape
    DG_col_edag(valid_edag_indices) = DG_col_edag_extracted(1:n_valid_modes); % Real parts
    DG_col_edag(n_edag + valid_edag_indices) = DG_col_edag_extracted(n_valid_modes+1:2*n_valid_modes); % Imaginary parts
    DG_col_edag(2*n_edag+1) = DG_col_edag_extracted(end-1); % Omega
    DG_col_edag(2*n_edag+2) = DG_col_edag_extracted(end); % Nu
    
    % Create identity vector: which Edaggershape mode does this QNshape column correspond to?
    % col_idx is the index in QNshape columns (1:DG_cols)
    % We need to check if this QNshape column corresponds to any Edaggershape mode
    evec1 = altzeros([A_size,1], nu);
    
    % Determine which QNshape nz=2 mode index col_idx corresponds to
    if col_idx <= n_vars_col  % Real part
        QN_mode_idx = col_idx;
        is_real = true;
    else  % Imaginary part
        QN_mode_idx = col_idx - n_vars_col;
        is_real = false;
    end
    
    % Find if this QNshape mode exists in Edaggershape
    % Compare with the modes in edag_to_QN_col_mapping
    edag_mode_match = find(edag_to_QN_col_mapping == QN_mode_idx);
    
    if ~isempty(edag_mode_match)
        % This QNshape column corresponds to an Edaggershape mode
        if is_real
            evec1(edag_mode_match) = 1; % Real part position in A matrix
        else
            evec1(n_edag + edag_mode_match) = 1; % Imaginary part position in A matrix
        end
    end
    % If no match, evec1 remains zero (this column is not in Edaggershape)

    % Compute I - A*DG for this column
    IADG_col = evec1 - A*DG_col_edag;

    % Use weights corresponding to Edaggershape structure
    weighted_IADGcol = abs(IADG_col).*weightsDG_edag;

    % Get the correct column weight 
    col_weight = weightsDG_col(col_idx);

    % Compute column sum and update maximum
    col_sum = sum(abs(weighted_IADGcol))/col_weight;
    max_col_sum = max(max_col_sum, col_sum);

    % Compute remainder contribution for this column
    if ~isempty(remainder_indices) && ~isempty(lambda_remainder)
        DG_col_remainder = DG_col(remainder_indices);
        evec_rem = altzeros([numel(DG_col_remainder),1], nu);
        match_rem = find(remainder_indices == col_idx);
        evec_rem(match_rem) = 1;
        Icol_remainder = abs(evec_rem - lambda_remainder .* remainder_weights .* DG_col_remainder);
        max_col_remainder = max(Icol_remainder) / col_weight;

        max_remainder_eigvec_cols = max(max_remainder_eigvec_cols, max_col_remainder);
    end
end

disp(['Max column sum: ', num2str(max_col_sum)]);

%% Compute omega and nu derivative columns

% Column for derivative w.r.t. omega
% From the DG structure: DG1 = Jr*x-Ji*y+omega*y  ->  d/d_omega = y
%                        DG2 = Ji*x+Jr*y-omega*x  ->  d/d_omega = -x
%                        DG3 = real(wT)*x-imag(wT)*y-1  ->  d/d_omega = 0
%                        DG4 = real(wT)*y+imag(wT)*x    ->  d/d_omega = 0

DG_col_omega = altzeros([DG_rows, 1], nu);
% Use nz=2 components from Q2Nshape structure - compute full column
nz2_positions = find(nz2_dg_row); % Positions of nz=2 components in Q2Nshape
x_nz2 = rv1large_Q2N(nz2_positions);  % real parts of nz=2 components
y_nz2 = iv1large_Q2N(nz2_positions);  % imaginary parts of nz=2 components

DG_col_omega(1:n_vars_row) = y_nz2;                    % DG1: d/d_omega = y
DG_col_omega(n_vars_row+1:2*n_vars_row) = -x_nz2;     % DG2: d/d_omega = -x
% DG3 and DG4 terms are zero for omega derivative

% Extract valid Edaggershape part for omega derivative and pad to full size
DG_col_omega_extracted = DG_col_omega(edag_row_indices);

% Create full-size omega column with zeros for missing modes
n_edag = sum(nz2_edag);
DG_col_omega_edag = altzeros([2*n_edag + 2, 1], nu);
n_valid_modes = sum(row_mapping_mask);
valid_edag_indices = find(row_mapping_mask); % Edaggershape mode indices that exist in Q2Nshape
DG_col_omega_edag(valid_edag_indices) = DG_col_omega_extracted(1:n_valid_modes); % Real parts
DG_col_omega_edag(n_edag + valid_edag_indices) = DG_col_omega_extracted(n_valid_modes+1:2*n_valid_modes); % Imaginary parts
DG_col_omega_edag(2*n_edag+1) = DG_col_omega_extracted(end-1); % Omega
DG_col_omega_edag(2*n_edag+2) = DG_col_omega_extracted(end); % Nu

% Compute I - A*DG_col_omega_edag
ID1 = eye(A_size);
IADG_col_omega = ID1(:,end-1)-A*DG_col_omega_edag;

% Compute weighted norm
weighted_IADGcol_omega = abs(IADG_col_omega).*weightsDG_edag;
col_sum_omega = sum(abs(weighted_IADGcol_omega))/etaOmega;
disp(['Omega column sum: ', num2str(col_sum_omega)]);

% Column for derivative w.r.t. nu

% Compute dudnu = -1/(2*nu^2)*classicforcingtensor
solshape.type='rec';
solshape.Nrec=N;

dudnu = -1/(2*nu^2)*classicforcingtensor;
dudnu = fulltosymmetrytensor_ext(setsizetensor(dudnu,N),solshape,symmetry);

% Set up reference solution
w0 = 1/(2*nu)*classicforcingtensor;
w1=fulltosymmetrytensor_ext(setsizetensor(w0,N),solshape,symmetry);
ref=w1(:);
ref=conjugatesymmetric_ext(ref,solshape,symmetry);

% Compute Jquad using FDFsymflextensor_ext_Hopf_trunc_nz2 with Q2Nshape to match DG structure
[~,~,~,Jquad]=FDFsymflextensor_ext_Hopf_trunc_nz2([dudnu(:);0],nu,classicforcingtensor,ref,solshape,symmetry,Q2Nshape);

N = sizeshape_Hopf(solshape);
M2N = max(sizeshape_Hopf(Q2Nshape), 2*N);

% Create the same index structure as in FDFsymflextensor_ext_Hopf_trunc_nz2
shapelarge.type = 'rec';
shapelarge.Nrec = M2N;
[nx_full, ny_full, nz_full, nt_full, ninshape_full, comp_full] = countersymtensor(M2N, shapelarge);

[tilden2_full, ~] = tildentensor(nx_full, ny_full, nz_full, nu);
% Compute symmetry indices
[symvar_mod, ~] = symmetryindicestensor_ext(nx_full, ny_full, nz_full, nt_full, comp_full, M2N, ninshape_full, symmetry);

% indices for the jacobian with Q2Nshape filtering
njac = shapetensor(nx_full, ny_full, nz_full, nt_full, Q2Nshape);
nz2_condition = (nz_full == 2);
jacsymvar = (symvar_mod & njac & nz2_condition);

DJnu_1 = diag(tilden2_full(jacsymvar));

% Combine to get full derivatives
DJr_nu = real(DJnu_1)+real(Jquad); % Derivative real(J) w.r.t. nu
DJi_nu = imag(DJnu_1)+imag(Jquad); % Derivative imag(J) w.r.t. nu

% Now compute the FULL nu column of DG
DG_col_nu = altzeros([DG_rows, 1], nu);

% Compute nu derivatives for the jacobian-sized components
DG1_nu_jac = DJr_nu*x_nz2 - DJi_nu*y_nz2;  % Real part derivative
DG2_nu_jac = DJi_nu*x_nz2 + DJr_nu*y_nz2;  % Imaginary part derivative

% Map to full Q2Nshape structure (fill the jacobian-computed parts)
DG_col_nu(1:n_vars_row) = DG1_nu_jac; % Real part derivatives
DG_col_nu(n_vars_row+1:n_vars_row+n_vars_row) = DG2_nu_jac; % Imaginary part derivatives
% DG3 and DG4 terms remain zero for nu derivative

% Extract valid Edaggershape part for nu derivative and pad to full size
DG_col_nu_extracted = DG_col_nu(edag_row_indices);

% Create fullsize nu column with zeros for missing modes
DG_col_nu_edag = altzeros([2*n_edag + 2, 1], nu);
valid_edag_indices = find(row_mapping_mask); % Edaggershape mode indices that exist in Q2Nshape
DG_col_nu_edag(valid_edag_indices) = DG_col_nu_extracted(1:n_valid_modes); % Real parts
DG_col_nu_edag(n_edag + valid_edag_indices) = DG_col_nu_extracted(n_valid_modes+1:2*n_valid_modes); % Imaginary parts
DG_col_nu_edag(2*n_edag+1) = DG_col_nu_extracted(end-1); % Omega
DG_col_nu_edag(2*n_edag+2) = DG_col_nu_extracted(end); % Nu

% Compute I - A*DG_col_nu_edag
IADG_col_nu = ID1(:,end)-A*DG_col_nu_edag;

% Compute weighted norm
weighted_IADGcol_nu = abs(IADG_col_nu).*weightsDG_edag;
col_sum_nu = sum(abs(weighted_IADGcol_nu))/etaNu;
disp(['Nu column sum: ', num2str(col_sum_nu)]);

% Update maximum with omega and nu contributions
max_col_sum = max([max_col_sum, col_sum_omega, col_sum_nu]);

%% Compute remainder contributions for omega and nu columns

disp(['Total remainder components processed: ', num2str(length(remainder_nx))]);
disp(['Eigvec remainder contribution (from loop): ', num2str(max_remainder_eigvec_cols)]);

% Extract remainder parts of omega and nu columns
DG_col_omega_rem= DG_col_omega(remainder_indices);
DG_col_nu_rem = DG_col_nu(remainder_indices);

% Create identity vectors for omega and nu (check if omega/nu rows are in remainder)
evec_rem_omega = altzeros([numel(DG_col_omega_rem),1], nu);
evec_rem_nu = altzeros([numel(DG_col_nu_rem),1], nu);

% Check if omega row (2*n_vars_row+1) is in remainder_indices
match_rem_omega = find(remainder_indices == 2*n_vars_row+1);
evec_rem_omega(match_rem_omega) = 1;

% Check if nu row (2*n_vars_row+2) is in remainder_indices  
match_rem_nu = find(remainder_indices == 2*n_vars_row+2);
evec_rem_nu(match_rem_nu) = 1;

% Compute weighted remainder contributions (including identity)
Icol_remainder_omega = abs(evec_rem_omega - lambda_remainder .* remainder_weights .* DG_col_omega_rem);
Icol_remainder_nu = abs(evec_rem_nu - lambda_remainder .* remainder_weights .* DG_col_nu_rem);

max_remainder_omega = max(Icol_remainder_omega) / etaOmega;
max_remainder_nu = max(Icol_remainder_nu) / etaNu;

disp(['Omega remainder contribution: ', num2str(max_remainder_omega)]);
disp(['Nu remainder contribution: ', num2str(max_remainder_nu)]);

% Total remainder contribution
max_remainder = max([max_remainder_eigvec_cols, max_remainder_omega, max_remainder_nu]);

disp(['Total remainder contribution: ', num2str(max_remainder)]);

% Return the maximum column sum
Z1finite = max(max_col_sum,max_remainder);
end