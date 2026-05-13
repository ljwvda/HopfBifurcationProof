function Z1finite = Z1finitebound_G_parallel_nz1(A,rv1large_Q2N,iv1large_Q2N,Omega,nu,wTrlarge_Q2N,wTilarge_Q2N,symmetry,Edaggershape,N,Etildeshape,eta,etaPhase)
% This is the nz=1 parallel implementation of Z1finitebound_G
% See Z1finitebound_G for detailed documentation of inputs and outputs.

% Set up extended shapes for finite bound computation
QNshape.type = 'otherplusrec';    % Column space
QNshape.other = Etildeshape;
QNshape.Nrec = N;

Q2Nshape.type = 'otherplusrec';   % Row space
Q2Nshape.other = Etildeshape;
Q2Nshape.Nrec = 2*N;

etaOmega = etaPhase(1);
etaNu = etaPhase(2);

%% Weights
% Use the finite part shape (Edaggershape) to compute weights
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(sizeshape_Hopf(Edaggershape),Edaggershape);
[symvar,~,~,grouporder,multiplicity,~] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,sizeshape_Hopf(Edaggershape),ninshape,symmetry);

%%

% Get nz=1 variables from the Q2Nshape space (for rows)
[nx_row,ny_row,nz_row,nt_row,ninshape_row,comp_row] = countersymtensor(sizeshape_Hopf(Q2Nshape),Q2Nshape);
[symvar_row,~,~,~,~,~] = symmetryindicestensor_ext(nx_row,ny_row,nz_row,nt_row,comp_row,sizeshape_Hopf(Q2Nshape),ninshape_row,symmetry);
nz1_row = (nz_row(symvar_row)==1);

% Get nz=1 variables from the QNshape space (for columns)
[nx_col,ny_col,nz_col,nt_col,ninshape_col,comp_col] = countersymtensor(sizeshape_Hopf(QNshape),QNshape);
[symvar_col,~,~,grouporder_col,multiplicity_col,~] = symmetryindicestensor_ext(nx_col,ny_col,nz_col,nt_col,comp_col,sizeshape_Hopf(QNshape),ninshape_col,symmetry);
nz1_col = (nz_col(symvar_col)==1);

% Compute weights for nz=1 column variables (QNshape)
orbits_col=grouporder_col./multiplicity_col; % orbit-stabilizer formula (from the DG col space)
weights_col=eta.^(abs(nx_col)+abs(ny_col)+abs(nz_col)+abs(nt_col)).*orbits_col;
weights_nz1_col = weights_col(symvar_col);
weights_nz1_col = weights_nz1_col(nz1_col); % Only nz=1 column variables

% For DG matrix:
% The DG matrix now has Q2Nshape rows and QNshape columns
% Rows: [x_nz1_row; y_nz1_row; omega; nu] where x_nz1_row, y_nz1_row are from Q2Nshape
% Cols: [x_nz1_col; y_nz1_col] where x_nz1_col, y_nz1_col are from QNshape
weightsDG_col = [weights_nz1_col; weights_nz1_col]; % Column weights
%%
% Create mapping to extract Edaggershape part from Q2Nshape structure
% We need to find which indices in Q2Nshape correspond to Edaggershape for nz=1 components

% Get Edaggershape indices for nz=1 components (this matches A matrix structure)
[nx_edag,ny_edag,nz_edag,nt_edag,ninshape_edag,comp_edag] = countersymtensor(sizeshape_Hopf(Edaggershape),Edaggershape);
[symvar_edag,~,~,grouporder_edag,multiplicity_edag,~] = symmetryindicestensor_ext(nx_edag,ny_edag,nz_edag,nt_edag,comp_edag,sizeshape_Hopf(Edaggershape),ninshape_edag,symmetry);
nz1_edag = (nz_edag(symvar_edag)==1);

n_vars_row = sum(nz1_row); % number of nz=1 row variables in Q2Nshape

% For each Edaggershape nz=1 mode, find its position in Q2Nshape

edag_to_Q2N_mapping = [];

% Get the (nx,ny,nz,nt) values for Edaggershape nz=1 components
nx_edag_nz1 = nx_edag(symvar_edag);
ny_edag_nz1 = ny_edag(symvar_edag);
nz_edag_nz1 = nz_edag(symvar_edag);
nt_edag_nz1 = nt_edag(symvar_edag);
comp_edag_nz1 = comp_edag(symvar_edag);

nx_edag_nz1 = nx_edag_nz1(nz1_edag);
ny_edag_nz1 = ny_edag_nz1(nz1_edag);
nz_edag_nz1 = nz_edag_nz1(nz1_edag);
nt_edag_nz1 = nt_edag_nz1(nz1_edag);
comp_edag_nz1 = comp_edag_nz1(nz1_edag);

% Get the (nx,ny,nz,nt,comp) values for Q2Nshape nz=1 components
nx_Q2N_nz1 = nx_row(symvar_row);
ny_Q2N_nz1 = ny_row(symvar_row);
nz_Q2N_nz1 = nz_row(symvar_row);
nt_Q2N_nz1 = nt_row(symvar_row);
comp_Q2N_nz1 = comp_row(symvar_row);

nx_Q2N_nz1 = nx_Q2N_nz1(nz1_row);
ny_Q2N_nz1 = ny_Q2N_nz1(nz1_row);
nz_Q2N_nz1 = nz_Q2N_nz1(nz1_row);
nt_Q2N_nz1 = nt_Q2N_nz1(nz1_row);
comp_Q2N_nz1 = comp_Q2N_nz1(nz1_row);

% Find matches: for each Edaggershape mode, find corresponding Q2Nshape position
for i = 1:length(nx_edag_nz1)
    match_idx = find(nx_Q2N_nz1 == nx_edag_nz1(i) & ...
        ny_Q2N_nz1 == ny_edag_nz1(i) & ...
        nz_Q2N_nz1 == nz_edag_nz1(i) & ...
        nt_Q2N_nz1 == nt_edag_nz1(i) & ...
        comp_Q2N_nz1 == comp_edag_nz1(i));
    if ~isempty(match_idx)
        edag_to_Q2N_mapping = [edag_to_Q2N_mapping; match_idx]; 
    else
        edag_to_Q2N_mapping = [edag_to_Q2N_mapping; -1];
    end
end

%% Create similar mapping for columns: QNshape to Edaggershape
% Find which QNshape nz=1 column indices correspond to Edaggershape modes
edag_to_QN_col_mapping = [];
n_vars_col = sum(nz1_col); % number of nz=1 column variables in QNshape

% Get the (nx,ny,nz,nt,comp) values for QNshape nz=1 column components
nx_QN_col_nz1 = nx_col(symvar_col);
ny_QN_col_nz1 = ny_col(symvar_col);
nz_QN_col_nz1 = nz_col(symvar_col);
nt_QN_col_nz1 = nt_col(symvar_col);
comp_QN_col_nz1 = comp_col(symvar_col);

nx_QN_col_nz1 = nx_QN_col_nz1(nz1_col);
ny_QN_col_nz1 = ny_QN_col_nz1(nz1_col);
nz_QN_col_nz1 = nz_QN_col_nz1(nz1_col);
nt_QN_col_nz1 = nt_QN_col_nz1(nz1_col);
comp_QN_col_nz1 = comp_QN_col_nz1(nz1_col);

% Find matches: for each Edaggershape mode, find corresponding QNshape column position
for i = 1:length(nx_edag_nz1)
    match_idx = find(nx_QN_col_nz1 == nx_edag_nz1(i) & ...
        ny_QN_col_nz1 == ny_edag_nz1(i) & ...
        nz_QN_col_nz1 == nz_edag_nz1(i) & ...
        nt_QN_col_nz1 == nt_edag_nz1(i) & ...
        comp_QN_col_nz1 == comp_edag_nz1(i));
    if ~isempty(match_idx)
        edag_to_QN_col_mapping = [edag_to_QN_col_mapping; match_idx]; 
    else
        edag_to_QN_col_mapping = [edag_to_QN_col_mapping; -1]; % Not found
    end
end

%%

% Create the index vector for [x_edag; y_edag; omega; nu]
% Filter to valid modes (exist in both Edaggershape and Q2Nshape)
row_mapping_mask = (edag_to_Q2N_mapping > 0);
valid_edag_to_Q2N_mapping = edag_to_Q2N_mapping(row_mapping_mask);

edag_row_indices = [valid_edag_to_Q2N_mapping; ... % x_edag (real parts)
    valid_edag_to_Q2N_mapping + n_vars_row; ... % y_edag (imaginary parts)
    2*n_vars_row+1; ... % omega
    2*n_vars_row+2]; % nu

% Weights for Edaggershape structure (matching A matrix)
orbits_edag=grouporder_edag./multiplicity_edag;
weights_edag=eta.^(abs(nx_edag)+abs(ny_edag)+abs(nz_edag)+abs(nt_edag)).*orbits_edag;
weights_nz1_edag = weights_edag(symvar_edag);
weights_nz1_edag = weights_nz1_edag(nz1_edag);
weightsDG_edag = [weights_nz1_edag; weights_nz1_edag; etaOmega; etaNu]; % [x; y; omega; nu] for Edaggershape
%%

A_size = size(A,1); % Size of A matrix (from Edaggershape)
DG_rows = 2*sum(nz1_row) + 2; % Size of DG rows (from Q2Nshape): 2*Q2N+2
DG_cols = 2*sum(nz1_col); % Size of DG columns (from QNshape): 2*QN

disp(['A matrix size: ', num2str(A_size), 'x', num2str(A_size)]);
disp(['DG matrix size: ', num2str(DG_rows), 'x', num2str(DG_cols)]);

% Pre-compute remainder indices and properties for use in the loop
all_indices = (1:DG_rows)';
remainder_indices = setdiff(all_indices, edag_row_indices); % Returns part of all_indices not in edag_row_indices

% Initialize remainder mode information
rem_nx = [];
rem_ny = [];
rem_nz = [];
rem_nt = [];
rem_comp = [];
rem_is_real = logical([]);  % true = real-part row, false = imaginary-part row

% Map remainder_indices back to the Q2Nshape structure to get (nx,ny,nz,nt,comp).
% DG row idx (1..n_vars_row) is the idx-th nz=1 mode in Q2Nshape symvar order,
% so directly index into the pre-computed nx_Q2N_nz1 / ny_Q2N_nz1 / ... arrays.
for i = 1:length(remainder_indices)
    idx = remainder_indices(i);
    if idx <= n_vars_row
        % Real-part row for the idx-th Q2Nshape nz=1 mode
        rem_nx = [rem_nx; nx_Q2N_nz1(idx)];
        rem_ny = [rem_ny; ny_Q2N_nz1(idx)];
        rem_nz = [rem_nz; nz_Q2N_nz1(idx)];
        rem_nt = [rem_nt; nt_Q2N_nz1(idx)];
        rem_comp = [rem_comp; comp_Q2N_nz1(idx)];
        rem_is_real = [rem_is_real; true];
    elseif idx <= 2*n_vars_row
        % Imaginary-part row for the (idx-n_vars_row)-th Q2Nshape nz=1 mode
        real_idx = idx - n_vars_row;
        rem_nx = [rem_nx; nx_Q2N_nz1(real_idx)];
        rem_ny = [rem_ny; ny_Q2N_nz1(real_idx)];
        rem_nz = [rem_nz; nz_Q2N_nz1(real_idx)];
        rem_nt = [rem_nt; nt_Q2N_nz1(real_idx)];
        rem_comp = [rem_comp; comp_Q2N_nz1(real_idx)];
        rem_is_real = [rem_is_real; false];
    end
    % Skip omega and nu terms (indices 2*n_vars_row+1 and 2*n_vars_row+2)
end

% Compute remainder lambda and weights
tilden2_remainder = rem_nx.^2 + rem_ny.^2 + rem_nz.^2;
lambda_remainder = 1./abs(nu*tilden2_remainder);
lambda_remainder(nu*tilden2_remainder==0)=0;
remainder_weights = eta.^(abs(rem_nx)+abs(rem_ny)+abs(rem_nz)+abs(rem_nt));

% Initialize arrays for parallel loop results
max_col_sum = altzeros([DG_cols, 1], nu); % Store column sums
max_remainder_eigvec_cols = altzeros([DG_cols, 1], nu); % Store remainder tracking for eigvec columns

% workaround to make parfor work with intlab
if exist('intval','file') && isintval(nu)
    workaroundintlabparfor(DG_cols);
end

disp('using parallelization for Z1finite bound computation (nz=1)');

% Loop over the DG columns - PARALLELIZED
parfor col_idx = 1:DG_cols

    if mod(col_idx, 1000) == 0
        disp(['Computing column ', num2str(col_idx), ' out of ', num2str(DG_cols)]);
    end

    % Compute single column of DG matrix with Q2Nshape structure
    DG_col = DGfunc_Hopf_nz1_column(rv1large_Q2N,iv1large_Q2N,Omega,nu,wTrlarge_Q2N,wTilarge_Q2N,N,symmetry,Q2Nshape,QNshape,col_idx);

    % Extract Edaggershape part from DG_col to match A matrix size
    DG_col_edag_extracted = DG_col(edag_row_indices); % Extract Edaggershape part
    
    % Create the full-size DG_col_edag vector with zeros for missing modes
    n_edag = sum(nz1_edag);  % Total number of Edaggershape nz=1 modes
    DG_col_edag = altzeros([2*n_edag + 2, 1], nu);
    
    % Fill in the valid parts (modes that exist in both Edaggershape and Q2Nshape)
    n_valid_modes = sum(row_mapping_mask);
    valid_edag_indices = find(row_mapping_mask); % Edaggershape mode indices that exist in Q2Nshape
    DG_col_edag(valid_edag_indices) = DG_col_edag_extracted(1:n_valid_modes); % Real parts
    DG_col_edag(n_edag + valid_edag_indices) = DG_col_edag_extracted(n_valid_modes+1:2*n_valid_modes); % Imaginary parts
    DG_col_edag(2*n_edag+1) = DG_col_edag_extracted(end-1); % Omega
    DG_col_edag(2*n_edag+2) = DG_col_edag_extracted(end); % Nu
    
    % Create identity vector: which Edaggershape mode does this QNshape column correspond to?
    evec1 = altzeros([A_size,1], nu);
    
    % Determine which QNshape nz=1 mode index col_idx corresponds to
    if col_idx <= n_vars_col  % Real part
        QN_mode_idx = col_idx;
        is_real = true;
    else  % Imaginary part
        QN_mode_idx = col_idx - n_vars_col;
        is_real = false;
    end
    
    % Find if this QNshape mode exists in Edaggershape
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

    IADG_col = evec1 - A*DG_col_edag;

    % Use weights corresponding to Edaggershape structure
    weighted_IADGcol = abs(IADG_col).*weightsDG_edag;

    col_weight = weightsDG_col(col_idx); % Use appropriate column weight (only x and y weights)

    % Compute column sum
    col_sum = sum(abs(weighted_IADGcol))/col_weight;
    % Workaround for intlab and parfor
    col_sum_scalar = col_sum(1,1);
    max_col_sum(col_idx) = col_sum_scalar;

    % Also compute remainder contribution for this column (using precomputed values)
    max_col_remainder = altzeros([1,1], nu); % Initialize
    if ~isempty(remainder_indices) && ~isempty(lambda_remainder)
        DG_col_remainder = DG_col(remainder_indices);
        evec_rem = altzeros([numel(DG_col_remainder),1], nu);
        % Find the remainder row whose mode and real/imaginary part matches this column.
        % col_idx indexes QNshape columns; remainder_indices index Q2Nshape rows.
        % These are different spaces, so compare by mode tuple (nx,ny,nz,nt,comp).
        if col_idx <= n_vars_col
            col_nx_val = nx_QN_col_nz1(col_idx);
            col_ny_val = ny_QN_col_nz1(col_idx);
            col_nz_val = nz_QN_col_nz1(col_idx);
            col_nt_val = nt_QN_col_nz1(col_idx);
            col_comp_val = comp_QN_col_nz1(col_idx);
            match_rem = find(rem_nx == col_nx_val & rem_ny == col_ny_val & ...
                rem_nz == col_nz_val & rem_nt == col_nt_val & ...
                rem_comp == col_comp_val & rem_is_real);
        else
            qn_idx = col_idx - n_vars_col;
            col_nx_val = nx_QN_col_nz1(qn_idx);
            col_ny_val = ny_QN_col_nz1(qn_idx);
            col_nz_val = nz_QN_col_nz1(qn_idx);
            col_nt_val = nt_QN_col_nz1(qn_idx);
            col_comp_val = comp_QN_col_nz1(qn_idx);
            match_rem = find(rem_nx == col_nx_val & rem_ny == col_ny_val & ...
                rem_nz == col_nz_val & rem_nt == col_nt_val & ...
                rem_comp == col_comp_val & ~rem_is_real);
        end
        evec_rem(match_rem) = 1;
        Icol_remainder = abs(evec_rem - lambda_remainder .* remainder_weights .* DG_col_remainder);
        max_col_remainder = max(Icol_remainder) / col_weight;
        max_col_remainder = max_col_remainder(1,1); % Workaround for intlab en parfor
    end

    max_remainder_eigvec_cols(col_idx) = max_col_remainder;


end

% After parallel loop, compute maximum values
max_col_sum = max(max_col_sum);
disp(['Max column sum: ', num2str(altsup(max_col_sum))]);
max_remainder_eigvec_cols = max(max_remainder_eigvec_cols);

%% Compute the last two columns of DG: derivatives w.r.t. omega and nu
% These are computed outside the for-loop as they don't depend on the column index

%% Column for derivative w.r.t. omega
% From the DG structure: DG1 = Jr*x-Ji*y+omega*y  ->  d/d_omega = y
%                        DG2 = Ji*x+Jr*y-omega*x  ->  d/d_omega = -x
%                        DG3 = real(wT)*x-imag(wT)*y-1  ->  d/d_omega = 0
%                        DG4 = real(wT)*y+imag(wT)*x    ->  d/d_omega = 0

DG_col_omega = altzeros([DG_rows, 1], nu);
% Use nz=1 components from Q2Nshape structure - compute FULL column
nz1_positions = find(nz1_row); % Positions of nz=1 components in Q2Nshape
x_nz1 = rv1large_Q2N(nz1_positions);  % real parts of nz=1 components
y_nz1 = iv1large_Q2N(nz1_positions);  % imaginary parts of nz=1 components

DG_col_omega(1:n_vars_row) = y_nz1;                    % DG1: d/d_omega = y
DG_col_omega(n_vars_row+1:2*n_vars_row) = -x_nz1;     % DG2: d/d_omega = -x
% DG3 and DG4 terms are zero for omega derivative

% Extract valid Edaggershape part for omega derivative and pad to full size
DG_col_omega_extracted = DG_col_omega(edag_row_indices);

% Create full-size omega column with zeros for missing modes
n_edag = sum(nz1_edag);
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
disp(['Omega column sum: ', num2str(altsup(col_sum_omega))]);


%% Column for derivative w.r.t. nu

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

% Compute Jquad using FDFsymflextensor_ext_Hopf_trunc_nz1 with Q2Nshape to match DG structure
[~,~,~,Jquad]=FDFsymflextensor_ext_Hopf_trunc_nz1([dudnu(:);0],nu,classicforcingtensor,ref,solshape,symmetry,Q2Nshape);

N = sizeshape_Hopf(solshape);
M2N = max(sizeshape_Hopf(Q2Nshape), 2*N);

% Create the same index structure as in FDFsymflextensor_ext_Hopf_trunc_nz1
shapelarge.type = 'rec';
shapelarge.Nrec = M2N;
[nx_full, ny_full, nz_full, nt_full, ninshape_full, comp_full] = countersymtensor(M2N, shapelarge);

[tilden2_full, ~] = tildentensor(nx_full, ny_full, nz_full, nu);
% Compute symmetry indices
[symvar_mod, ~] = symmetryindicestensor_ext(nx_full, ny_full, nz_full, nt_full, comp_full, M2N, ninshape_full, symmetry);

% indices for the jacobian with Q2Nshape filtering
njac = shapetensor(nx_full, ny_full, nz_full, nt_full, Q2Nshape);
nz1_condition = (nz_full == 1);
jacsymvar = (symvar_mod & njac & nz1_condition);

DJnu_1 = diag(tilden2_full(jacsymvar));

% Combine to get full derivatives
DJr_nu = real(DJnu_1)+real(Jquad); % Derivative real(J) w.r.t. nu
DJi_nu = imag(DJnu_1)+imag(Jquad); % Derivative imag(J) w.r.t. nu

% Now compute the full nu column of DG
DG_col_nu = altzeros([DG_rows, 1], nu);

% Compute nu derivatives for the jacobian-sized components
DG1_nu_jac = DJr_nu*x_nz1 - DJi_nu*y_nz1;  % Real part derivative
DG2_nu_jac = DJi_nu*x_nz1 + DJr_nu*y_nz1;  % Imaginary part derivative

% Map to FULL Q2Nshape structure
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

disp(['Nu column sum: ', num2str(altsup(col_sum_nu))]);


%% Update maximum with omega and nu contributions
max_col_sum = max([max_col_sum, col_sum_omega, col_sum_nu]);

%% Compute remainder contributions (remainder indices and modes info computed before loop)

disp(['Total remainder components processed: ', num2str(length(rem_nx))]);

disp(['eigvec remainder contribution (from loop): ', num2str(altsup(max_remainder_eigvec_cols))]);


% Now compute remainder contributions for omega and nu derivative columns
% Extract remainder parts of omega and nu columns
DG_col_omega_remainder = DG_col_omega(remainder_indices);
DG_col_nu_remainder = DG_col_nu(remainder_indices);

% Create identity vectors for omega and nu (check if omega/nu rows are in remainder)
evec_rem_omega = altzeros([numel(DG_col_omega_remainder),1], nu);
evec_rem_nu = altzeros([numel(DG_col_nu_remainder),1], nu);

% Check if omega row (2*n_vars_row+1) is in remainder_indices
match_rem_omega = find(remainder_indices == 2*n_vars_row+1);
evec_rem_omega(match_rem_omega) = 1;

% Check if nu row (2*n_vars_row+2) is in remainder_indices  
match_rem_nu = find(remainder_indices == 2*n_vars_row+2);
evec_rem_nu(match_rem_nu) = 1;

% Compute weighted remainder contributions (including identity)
Icol_remainder_omega = abs(evec_rem_omega - lambda_remainder .* remainder_weights .* DG_col_omega_remainder);
Icol_remainder_nu = abs(evec_rem_nu - lambda_remainder .* remainder_weights .* DG_col_nu_remainder);

max_remainder_omega = max(Icol_remainder_omega) / etaOmega;
max_remainder_nu = max(Icol_remainder_nu) / etaNu;

disp(['Omega remainder contribution: ', num2str(altsup(max_remainder_omega))]);
disp(['Nu remainder contribution: ', num2str(altsup(max_remainder_nu))]);

% Total remainder contribution
max_remainder = max([max_remainder_eigvec_cols, max_remainder_omega, max_remainder_nu]);

disp(['Total remainder contribution: ', num2str(altsup(max_remainder))]);


% Return the maximum column sum
Z1finite = max(max_col_sum,max_remainder);

end
