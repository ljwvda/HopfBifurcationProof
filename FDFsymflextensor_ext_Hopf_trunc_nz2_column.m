function [J_col] = FDFsymflextensor_ext_Hopf_trunc_nz2_column(x,nu,~,~,solshape,symmetry,Jshape,col_idx)
% Computes a single column of the Jacobian
% 
% This is an optimized version of FDFsymflextensor_ext_Hopf_trunc_nz2 that computes only
% the specified column of the Jacobian matrix without building the full
% matrix.
%
% Inputs are the same as FDFsymflextensor_ext_Hopf_trunc_nz2, plus:
%   col_idx - Column index to compute
%
% Output:
%   J_col - Single column of the Jacobian matrix

%%%%%%%%%%%%%%%%%% INDEXING %%%%%%%%%%%%%%%

if ~exist('Jshape','var')
    Jshape=solshape;
end

%%%%%%%%%%%%%%%%%% INDEXING %%%%%%%%%%%%%%%

if ~exist('Jshape','var')
    Jshape=solshape;
end
N=sizeshape_Hopf(solshape);
M2N=max(sizeshape_Hopf(Jshape),2*N);

% biggest index set needed
shapelarge.type='rec';
shapelarge.Nrec=M2N;
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(M2N,shapelarge);

% Compute symmetry indices with original nz for compatibility
[symvar_orig,~] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M2N,ninshape,symmetry);
[tilden2,tilden2reci] = tildentensor(nx,ny,nz,nu);

% indices for the solution
nsol=shapetensor(nx,ny,nz,nt,solshape);
solsymvar=(symvar_orig & nsol);
[solnx,solny,solnz] = countersymtensor(N,solshape);
soltilden2reci=setsizetensor(tilden2reci,N);

% Apply ellHopf condition for Jacobian computation
if strcmp(Jshape.type, 'ellHopf')% Specific for Hopf bifurcation problem with nz=2
    nz = 2*ones(size(nz));
end

% Recompute symmetry indices with modified nz and modified symmetry 
symmetry_idx = symmetry;
if (strcmp(Jshape.type, 'ellHopf') && symmetry == 1)
    symmetry_idx = 11;  
elseif (strcmp(Jshape.type, 'ellHopf') && symmetry == 25)
    symmetry_idx = 26;  
end

[symvar,symindex] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M2N,ninshape,symmetry_idx);

% indices for the jacobian 
njac = shapetensor(nx,ny,nz,nt,Jshape);
nz2_condition = (nz == 2);
jacsymvar = (symvar & njac & nz2_condition);

MJ = length(find(jacsymvar));

% Validate column index
if col_idx < 1 || col_idx > MJ
    error('Column index %d is out of bounds. Matrix has %d columns.', col_idx, MJ);
end

% extract solution coefficients
wsol=x(1:end-1);

% Note: wjac computation is removed as it's not needed for column-wise computation

%%%%%%%%%% Computation of linear and nonlinear part %%%%%%%%%%

% From wsol, build a tensor of the right size and convert it to symmetrized variables
w=symmetrytofulltensor_ext(wsol,solshape,symmetry);
w=setsizetensor(w,M2N); 
% Extract jacobian part (though we don't use wjac directly in column computation)
% wjac=w(jacsymvar); 
w=setsizetensor(w,N); % Resize back to N for derivative computation

% derivative operators
Dw=altzeros([2*N+1,3,3],nu);
Dw(:,:,:,:,:,1)=P(solnx,w);
Dw(:,:,:,:,:,2)=P(solny,w);
Dw(:,:,:,:,:,3)=P(solnz,w);

% omega cross w  
Mw=altzeros([2*N+1,3],nu);
matrix=[0 -3 2; 3 0 -1; -2 1 0]; 
for k=1:3
   vec=matrix(k,:);
   coef=find(vec~=0);
   Mw(:,:,:,:,k)=sign(vec(coef(1)))*Dw(:,:,:,:,coef(1),abs(vec(coef(1))))+...
                   sign(vec(coef(2)))*Dw(:,:,:,:,coef(2),abs(vec(coef(2))));
end

Mw=P(Mw,1i*soltilden2reci);

% the derivative of M*omega
DMw=altzeros([2*N+1,3,3],nu);
DMw(:,:,:,:,:,1)=P(solnx,Mw);
DMw(:,:,:,:,:,2)=P(solny,Mw);
DMw(:,:,:,:,:,3)=P(solnz,Mw);

%%%%%%%%%% JACOBIAN COLUMN COMPUTATION %%%%%%%%%%

% Compute the modified symmetry index 
idx_symmetry = symmetry;
if (strcmp(Jshape.type, 'ellHopf') && symmetry == 1)
    idx_symmetry = 11;  
elseif (strcmp(Jshape.type, 'ellHopf') && symmetry == 25)
    idx_symmetry = 26;  
end

% derivative of the nonlinear part in omega - compute only the specified column
J_col = Dnonlinear_nz2_column(w,Dw,Mw,DMw,N,Jshape,idx_symmetry,col_idx);

% Add diagonal term for the Laplacian (only affects the col_idx-th element)
auxJ = nu*tilden2(jacsymvar);
auxJ_vec = auxJ(:);
J_col(col_idx) = J_col(col_idx) + auxJ_vec(col_idx);

end  % of main function 

function J_col=Dnonlinear_nz2_column(w,Dw,Mw,DMw,N,Jshape,symmetry,col_idx)
% determines a single column of the Jacobian of the nonlinear part 
% This version computes ONLY the requested column without building the full Jacobian matrix

  % initialize indices - use same approach as full function but filter for nz=2
  M=sizeshape_Hopf(Jshape); 
  MN=max(M,N);
  [nx,ny,nz,nt,ninshape,comp] = countersymtensor(MN,Jshape);
  if strcmp(Jshape.type, 'ellHopf') % Specific for Hopf bifurcation problem with nz=2
      nz = 2*ones(size(nz));
  end
  [symvar,symindex,symfactor] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,MN,ninshape,symmetry);
  
  % Filter for nz == 2 only
  nz2_condition = (nz == 2);
  symvar_nz2 = symvar & nz2_condition;
  
  [~,tilden2reci] = tildentensor(nx,ny,nz,Dw(1));
  JM=length(find(symvar_nz2));
  J_col=altzeros([JM,1], Dw(1));  % Initialize column vector with interval compatibility
  [nsymvar,NN]=makeindices(N,nx,ny,nz,nt,symvar); 
  
  % Create mapping from full symvar indices to nz2 symvar indices
  full_symvar_indices = find(symvar);
  nz2_symvar_indices = find(symvar_nz2);
  [~, nz2_row_mapping] = ismember(nz2_symvar_indices, full_symvar_indices);
  
  % Create the inverse mapping: from full symvar positions to nz2 positions
  full_to_nz2_mapping = zeros(length(full_symvar_indices), 1);  
  for i = 1:length(nz2_row_mapping)
    if nz2_row_mapping(i) > 0
      full_to_nz2_mapping(nz2_row_mapping(i)) = i;
    end
  end
  
  compsymvar=comp(symvar);
  
  % Loop over all jjj that have nz=2 (same as full function)
  for jjj = find(symindex & nz2_condition)'
    
    % Find which column this jjj contributes to
    j = symindex(jjj);  % This gives the reduced symmetry index
    [~, col_index] = ismember(j, nz2_row_mapping);
    
    % Only process if this jjj contributes to our requested column
    if col_index ~= col_idx
      continue;
    end
    
    % identify the nonzero components of M_n^{l,m}
    matrix=[0 -3 2; 3 0 -1; -2 1 0]; 
    column=matrix(:,comp(jjj)); 
    lvals=find(column~=0); %nonvanishing components l in M_n^{l,m}
    nvals=column(lvals(:));
    nnreci=1i*[nx(jjj);ny(jjj);nz(jjj)]*tilden2reci(jjj);  
    Mvals=sign(nvals).*nnreci(abs(nvals));

    % shift of the indices (coming from derivative of convolution)
    nsx=nsymvar(:,1)-nx(jjj);
    nsy=nsymvar(:,2)-ny(jjj);
    nsz=nsymvar(:,3)-nz(jjj);
    nst=nsymvar(:,4)-nt(jjj);
    
    % only the shifted indices within the numerical solution rectangle 
    % will lead to a nonvanishing contribution, AND only nz=2 elements
    jj=(abs(nsx)<=N(1) & abs(nsy)<=N(2) & abs(nsz)<=N(3) & abs(nst)<=N(4)) & (nsymvar(:,3) == 2);
    
    if ~any(jj)
        continue;
    end
    
    % convert tensor indices into linear indices for nx,ny,nz dimensions of tensor
    Sall1=1+nsx(jj)+N(1)+NN(1)*(nsy(jj)+N(2))+NN(2)*(nsz(jj)+N(3))+NN(3)*(nst(jj)+N(4));
    % next, including the nt dimension
    Sall2=Sall1+NN(4)*(compsymvar(jj)-1); 
    % including the component dimension, but selecting the component corresponding to jjj
    Sall3=Sall2+NN(5)*(comp(jjj)-1);

    % Compute the four terms in DPsi
    
    % term sum_{p=1}^3 n_p (Mw)^p_{k-n}
    DPsi1=Mw(Sall1)*nx(jjj)+Mw(Sall1+NN(4))*ny(jjj)+Mw(Sall1+2*NN(4))*nz(jjj);
    % the factor delta_{l,m}
    DPsi1=DPsi1.*(compsymvar(jj)==comp(jjj));
    
    % term (D_m(Mw)^l)_{k-n}
    DPsi2=DMw(Sall3);

    % term sum_{p=1}^3 M_n^{p,m} (D_p w^l)_{k-n}
    DPsi3=Dw(Sall2+NN(5)*(lvals(1)-1))*Mvals(1)+Dw(Sall2+NN(5)*(lvals(2)-1))*Mvals(2);          

    % term sum_{p=1}^3 n_p w^p_{k-n}
    DPsi4=w(Sall1)*nx(jjj)+w(Sall1+NN(4))*ny(jjj)+w(Sall1+2*NN(4))*nz(jjj);
    % the factor M_n^{l,m}
    Mnlm_expanded = Mvals(1)*(compsymvar(jj)==lvals(1))+Mvals(2)*(compsymvar(jj)==lvals(2));
    DPsi4=DPsi4.*Mnlm_expanded;
      
    % Map to nz2-only indices and add contributions to the column
    % find(jj) gives positions in the full symvar list (via nsymvar)
    % We need to map these to nz2-only row indices
    jj_positions = find(jj);
    row_indices = full_to_nz2_mapping(jj_positions);
    
    valid_mask = row_indices > 0;
    if any(valid_mask)
      valid_rows = row_indices(valid_mask);
      J_col(valid_rows) = J_col(valid_rows) + symfactor(jjj)*(DPsi1(valid_mask)-DPsi2(valid_mask)+DPsi3(valid_mask)-DPsi4(valid_mask));
    end
    
  end  % End of jjj loop

  % include the factor i
  J_col=1i*J_col;
end

%%%%%%%%%%%%%%%%%%%%%% Some auxiliary functions %%%%%%%%%%%%%%%%%%

function product=P(a,b)
% elementwise product of tensors for compatibility with intlab
    product=reshape(a(:).*b(:),size(a));    
end
