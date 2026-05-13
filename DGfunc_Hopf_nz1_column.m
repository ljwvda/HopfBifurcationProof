function DG_col = DGfunc_Hopf_nz1_column(x,y,omega,nu,wTr,wTi,N,symmetry,Q2Nshape,QNshape,col_idx)
%% Computes a single column of the DG matrix for the Hopf bifurcation system (nz=1)
% This is the column-wise version of DGfunc_Hopf_nz1.m
% Modified to use Q2Nshape (for rows) and QNshape (for columns) instead of Edaggershape
%
% Inputs:
%   x, y - eigenvectors (not yet filtered for nz=1)
%   omega - frequency
%   nu - viscosity 
%   wTr, wTi - real and imaginary parts of wT (not yet filtered for nz=1)
%   N - truncation parameters [N1, N2, N3, N4]
%   symmetry - symmetry parameter
%   Q2Nshape - shape structure for DG matrix rows (size 2*Q2N+2)
%   QNshape - shape structure for DG matrix columns (size 2*QN)
%   col_idx - column index to compute 
%
% Output:
%   DG_col - single column of the DG matrix (size 2*Q2N+2 by 1)

solshape.type='rec';
solshape.Nrec=N;

%% Setup reference solution
forcing = classicforcingtensor;
if exist('intval','file') && isintval(nu)
    w0 = 1/(2*nu)*forcing; 
    w0 = intval(w0);
else
    w0 = 1/(2*nu)*forcing; 
end
w1=fulltosymmetrytensor_ext(setsizetensor(w0,N),solshape,symmetry);
ref=w1(:); % reference solution
ref=conjugatesymmetric_ext(ref,solshape,symmetry);

%% Get matrix dimensions using Q2Nshape for rows and QNshape for columns
% For column variables (QNshape): biggest index set needed for column space
[nx_col,ny_col,nz_col,nt_col,ninshape_col,comp_col] = countersymtensor(sizeshape_Hopf(QNshape),QNshape);
[symvar_col,~] = symmetryindicestensor_ext(nx_col,ny_col,nz_col,nt_col,comp_col,sizeshape_Hopf(QNshape),ninshape_col,symmetry);

% For row variables (Q2Nshape): biggest index set needed for row space
[nx_row,ny_row,nz_row,nt_row,ninshape_row,comp_row] = countersymtensor(sizeshape_Hopf(Q2Nshape),Q2Nshape);
[symvar_row,~] = symmetryindicestensor_ext(nx_row,ny_row,nz_row,nt_row,comp_row,sizeshape_Hopf(Q2Nshape),ninshape_row,symmetry);

% Filter for nz=1 components
nz1_col = (nz_col(symvar_col)==1);
nz1_row = (nz_row(symvar_row)==1);

% The input vectors x, y, wTr, wTi have size corresponding to Q2Nshape
x_row = x(nz1_row);
y_row = y(nz1_row);
wTr_row = wTr(nz1_row);
wTi_row = wTi(nz1_row);

n_vars_col = sum(nz1_col); % number of nz=1 column variables (QNshape)
n_vars_row = sum(nz1_row); % number of nz=1 row variables (Q2Nshape)

%% Create mapping from QNshape column indices to Q2Nshape row indices
% Get (nx,ny,nz,nt) for QNshape nz=1 modes
nx_col_nz1 = nx_col(symvar_col);
ny_col_nz1 = ny_col(symvar_col);
nz_col_nz1 = nz_col(symvar_col);
nt_col_nz1 = nt_col(symvar_col);
comp_col_nz1 = comp_col(symvar_col);

nx_col_nz1 = nx_col_nz1(nz1_col);
ny_col_nz1 = ny_col_nz1(nz1_col);
nz_col_nz1 = nz_col_nz1(nz1_col);
nt_col_nz1 = nt_col_nz1(nz1_col);
comp_col_nz1 = comp_col_nz1(nz1_col);

% Get (nx,ny,nz,nt) for Q2Nshape nz=1 modes
nx_row_nz1 = nx_row(symvar_row);
ny_row_nz1 = ny_row(symvar_row);
nz_row_nz1 = nz_row(symvar_row);
nt_row_nz1 = nt_row(symvar_row);
comp_row_nz1 = comp_row(symvar_row);

nx_row_nz1 = nx_row_nz1(nz1_row);
ny_row_nz1 = ny_row_nz1(nz1_row);
nz_row_nz1 = nz_row_nz1(nz1_row);
nt_row_nz1 = nt_row_nz1(nz1_row);
comp_row_nz1 = comp_row_nz1(nz1_row);

if isequal(Q2Nshape, QNshape) && n_vars_col == n_vars_row
    % Verify that the nz=1 modes are ordered identically
    if isequal(nx_col_nz1, nx_row_nz1) && isequal(ny_col_nz1, ny_row_nz1) && ...
       isequal(nz_col_nz1, nz_row_nz1) && isequal(nt_col_nz1, nt_row_nz1)
        % Use identity mapping for efficiency and correctness
        col_to_row_map = (1:n_vars_col)';
    else
        error('Q2Nshape and QNshape are identical but nz=1 modes differ in ordering');
    end
else
    col_to_row_map = zeros(n_vars_col, 1);
    for i = 1:n_vars_col
                match_idx = find(nx_row_nz1 == nx_col_nz1(i) & ...
                         ny_row_nz1 == ny_col_nz1(i) & ...
                         nz_row_nz1 == nz_col_nz1(i) & ...
                         nt_row_nz1 == nt_col_nz1(i) & ...
                         comp_row_nz1 == comp_col_nz1(i));
        if ~isempty(match_idx)
            col_to_row_map(i) = match_idx;
        else
            error('Column mode %d not found in row space', i);
        end
    end
end

% Initialize output column with proper size (2*Q2N+2 rows)
DG_col = altzeros([2*n_vars_row + 2, 1], nu);

%% Compute derivatives based on col_idx
if col_idx <= n_vars_col
    % Derivative w.r.t. x(col_idx) - real part of QNshape variable
    var_idx_col = col_idx;  % Index in QNshape column space
    var_idx_row = col_to_row_map(var_idx_col);  % Corresponding index in Q2Nshape row space
    
    % Compute Jacobian column for the correct variable (using Q2Nshape indexing)
    J_col = FDFsymflextensor_ext_Hopf_trunc_nz1_column([w1(:);0],nu,forcing,ref,solshape,symmetry,Q2Nshape,var_idx_row);
    Jr_col = real(J_col);
    Ji_col = imag(J_col);
    
    % DG1 = Jr*x-Ji*y+omega*y  ->  d/dx_k = Jr(:,k)
    DG_col(1:n_vars_row) = Jr_col;
    
    % DG2 = Ji*x+Jr*y-omega*x  ->  d/dx_k = Ji(:,k) - omega*delta_k
    DG_col(n_vars_row+1:2*n_vars_row) = Ji_col;
    % Add omega term at the correct row position (var_idx_row, not var_idx_col)
    DG_col(n_vars_row+var_idx_row) = DG_col(n_vars_row+var_idx_row) - omega;
    
    % DG3 = real(wT)*x-imag(wT)*y-1  ->  d/dx_k = wTr(k)
    DG_col(2*n_vars_row+1) = wTr_row(var_idx_row);
    
    % DG4 = real(wT)*y+imag(wT)*x  ->  d/dx_k = wTi(k)
    DG_col(2*n_vars_row+2) = wTi_row(var_idx_row);
    
elseif col_idx <= 2*n_vars_col
    % Derivative w.r.t. y(col_idx - n_vars_col) - imaginary part of QNshape variable
    var_idx_col = col_idx - n_vars_col;  % Index in QNshape column space
    var_idx_row = col_to_row_map(var_idx_col);  % Corresponding index in Q2Nshape row space
    
    % Compute Jacobian column  (using Q2Nshape indexing)
    J_col = FDFsymflextensor_ext_Hopf_trunc_nz1_column([w1(:);0],nu,forcing,ref,solshape,symmetry,Q2Nshape,var_idx_row);
    Jr_col = real(J_col);
    Ji_col = imag(J_col);
    
    % DG1 = Jr*x-Ji*y+omega*y  ->  d/dy_k = -Ji(:,k) + omega*delta_k  
    DG_col(1:n_vars_row) = -Ji_col;
    % Add omega term at the correct row position (var_idx_row, not var_idx_col)
    DG_col(var_idx_row) = DG_col(var_idx_row) + omega;
    
    % DG2 = Ji*x+Jr*y-omega*x  ->  d/dy_k = Jr(:,k)
    DG_col(n_vars_row+1:2*n_vars_row) = Jr_col;
    
    % DG3 = real(wT)*x-imag(wT)*y-1  ->  d/dy_k = -wTi(k)
    DG_col(2*n_vars_row+1) = -wTi_row(var_idx_row);
    
    % DG4 = real(wT)*y+imag(wT)*x  ->  d/dy_k = wTr(k)
    DG_col(2*n_vars_row+2) = wTr_row(var_idx_row);
    
else
    error('Column index out of range');
end

end
