function [F,J,DFnu,Jquad,Fext,Fphase,tensors,F_nz2] = ...
             FDFsymflextensor_ext_Hopf_trunc_nz2(x,nu,forcing,ref,solshape,symmetry,Jshape)
% Optimized version that only computes Jacobian for nz==2 elements
% x is assumed to be in symmetrized variables form (of shape solshape)
% g is assumed to be in (perhaps small) tensor form
% nu is the parameter
% ref is the reference solution (symmetrized) fixing the time shift
% Jshape determines the indices included in the Jacobian (default Jshape=solshape)
%
% F is the residue (n symmetrized variables) of length sizeshape(solshape)+1
% J is the Jacobian (w.r.t. x) of size sizeshape(Jshape)+1 
% DFnu is the derivative of F w.r.t. nu
% Jquad is the "quadratic part" (useful for bifurcation problems)
% i.e. the part of the Jacobian that is non-constant (depends on x)
% Fext is the extended residue (all nonzero coefficients) 
% as a tensor of size 2N (no phase eqn), where N=sizeshape(solshape)
% Fphase is the phase condition
% tensors contains additional output arguments: Dw,Mw,DMw as tensors of size N
% F_nz2 is the residue filtered for only nz=2 components (without phase condition)

% Modified version of FDFsymflextensor_ext for Hopf bifurcation, excluding
% the phase condition in the Jacobian matrix. Only computes nz==2 elements.

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
[symvar_orig,symindex_orig] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M2N,ninshape,symmetry);
[tilden2,tilden2reci] = tildentensor(nx,ny,nz,nu);

% indices for the solution
nsol=shapetensor(nx,ny,nz,nt,solshape);
solsymvar=(symvar_orig & nsol);
symindexsol=symindex_orig(solsymvar);
[solnx,solny,solnz] = countersymtensor(N,solshape);
soltilden2reci=setsizetensor(tilden2reci,N);

% Apply ellHopf condition for Jacobian computation
if strcmp(Jshape.type, 'ellHopf')% Specific for Hopf bifurcation problem with nz=2
    nz = 2*ones(size(nz));
end

% Recompute symmetry indices with modified nz 
% For ellHopf: use symmetry=11 to include components 1 and 2 for nz=2, but only for Jacobian
idx_symmetry = symmetry;
if (strcmp(Jshape.type, 'ellHopf') & symmetry == 1)
    idx_symmetry = 11;  % This ensures we include components 1 and 2 for nz=2
elseif (strcmp(Jshape.type, 'ellHopf') & symmetry == 25)
    idx_symmetry = 26;  % This ensures we include components 1 and 2 for nz=2
end
[symvar,symindex] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M2N,ninshape,idx_symmetry);

% indices for the jacobian - use the modified symmetry indices
njac = shapetensor(nx,ny,nz,nt,Jshape);
nz2_condition = (nz == 2);
jacsymvar = (symvar & njac & nz2_condition);
symindexjac = symindex(jacsymvar);

% Recompute tilden2 with modified nz (for both residual and Jacobian)
[tilden2,tilden2reci] = tildentensor(nx,ny,nz,nu); 

% symvar_nz2 = symvar & nz2_condition;
 
wsol=x(1:end-1);
w=symmetrytofulltensor_ext(wsol,solshape,symmetry);
w=setsizetensor(w,M2N); 
wjac=w(jacsymvar);
w=setsizetensor(w,N);

Omega=x(end);

% extend the reference for the phase condition
reflarge=altzeros([length(find(symvar_orig)),1],nu);
reflarge(symindexsol)=ref;
% refjac=reflarge(symindexjac); 

%%%%%%%%%%%%%%%%%% FUNCTION F %%%%%%%%%%%%%%%%%%%

% derivative of omega
Dw=altzeros([2*N+1,3,3],nu);
Dw(:,:,:,:,:,1)=P(solnx,w); 
Dw(:,:,:,:,:,2)=P(solny,w);
Dw(:,:,:,:,:,3)=P(solnz,w);

% determine M*omega
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

if nargout>=7
    tensors.Dw=Dw;
    tensors.Mw=Mw;
    tensors.DMw=DMw;
end

% calculate F
forcingsol=fulltosymmetrytensor_ext(setsizetensor(forcing,N),solshape,symmetry);
F1 = (nu*tilden2(solsymvar)).*wsol - forcingsol; 

F1large2N=setsizetensor(symmetrytofulltensor_ext(F1,solshape,symmetry),2*N);
F2large2N=tensorproductext(Mw,Dw); 
F3large2N=tensorproductext(w,DMw); 
Fext=F1large2N+1i*(F2large2N-F3large2N); % this has tensor data structure

Fsol=fulltosymmetrytensor_ext(setsizetensor(Fext,N),solshape,symmetry);

% phase condition
Fphase=zeros(size(1i*ref'*(wsol))); 

F=[Fsol;Fphase];

% Compute F filtered for nz=2 components only (without phase condition)
% Use the exact same approach as the Jacobian - extract from full tensor first, then symmetrize
Fext_truncated = setsizetensor(Fext, M2N);
F_nz2 = Fext_truncated(jacsymvar); 

if nargout==1 % no Jacobian computed if not needed
    return
end


%%%%%%%%%% JACOBIAN %%%%%%%%%%%% 

% derivative of the nonlinear part in omega - for nz==2 only
Jnonlinear=Dnonlinear_nz2(w,Dw,Mw,DMw,N,Jshape,idx_symmetry);

MJ=length(find(jacsymvar));
J=altzeros([MJ,MJ],nu);
J(1:MJ,1:MJ)=Jnonlinear; 

diagMJ=(speye(MJ)==1);

% derivative of the time derivative term w.r.t. omega
J(diagMJ)=J(diagMJ);
% Jquad is an extra output containing the derivative of all quadratic terms
Jquad=J; 
% derivative of the Laplacian term
auxJ = nu*tilden2(jacsymvar);
J(diagMJ)=J(diagMJ)+auxJ;

% No phase condition derivatives needed for this truncated version

% DFnu is an extra output containing the derivative w.r.t. the parameter nu
DFnu=tilden2(jacsymvar).*wjac;

end  % of main function 

function JJ=Dnonlinear_nz2(w,Dw,Mw,DMw,N,Jshape,symmetry)
% determines the Jacobian of the nonlinear part - FOR NZ==2 ONLY

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
  JJ=altzeros([JM,JM],Dw(1));
  [nsymvar,NN]=makeindices(N,nx,ny,nz,nt,symvar); 
  
  % Create mapping from full symvar indices to nz2 symvar indices
  full_symvar_indices = find(symvar);
  nz2_symvar_indices = find(symvar_nz2);
  [~, nz2_row_mapping] = ismember(nz2_symvar_indices, full_symvar_indices);
  
  compsymvar=comp(symvar);
  compsymvar_nz2=comp(symvar_nz2);
  compind{1}=(compsymvar_nz2==1);
  compind{2}=(compsymvar_nz2==2);
  compind{3}=(compsymvar_nz2==3);

  % for loop to fill the columns of the matrix - only iterate over nz == 2 elements
  for jjj=find(symindex & nz2_condition)'      
    % for jjj=find(symindex)'      
    % derivative w.r.t. jjj-th symmetry (not necessarily reduced) variable 
    
    % identify the nonzero components of M_n^{l,m}
    % needed for the terms DPsi2 and DPsi3
    matrix=[0 -3 2; 3 0 -1; -2 1 0]; 
    column=matrix(:,comp(jjj)); 
    lvals=find(column~=0); %nonvanishing components l in M_n^{l,m}
    nvals=column(lvals(:));
    nnreci=1i*[nx(jjj);ny(jjj);nz(jjj)]*tilden2reci(jjj);  
    Mvals=sign(nvals).*nnreci(abs(nvals));
    Mnlm=Mvals(1)*(compind{lvals(1)})+Mvals(2)*(compind{lvals(2)});
 
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
    % as this is the one needed below in DPsi3
    Sall3=Sall2+NN(5)*(comp(jjj)-1);

    % Compute the four terms in DPsi
    
    % term sum_{p=1}^3 n_p (Mw)^p_{k-n}
    % note that sum_{p=1}^3 n_p(Mw)^p_{k-n} = sum_{p=1}^3 k_p(Mw)^p_{k-n}
    % since Mw is divergence free for any w
    DPsi1=Mw(Sall1)*nx(jjj)+Mw(Sall1+NN(4))*ny(jjj)+Mw(Sall1+2*NN(4))*nz(jjj);
    % the factor delta_{l,m}
    DPsi1=DPsi1.*(compsymvar(jj)==comp(jjj));
    
    % term (D_m(Mw)^l)_{k-n}
    DPsi2=DMw(Sall3);

    % term sum_{p=1}^3 M_n^{p,m} (D_p w^l)_{k-n}
    DPsi3=Dw(Sall2+NN(5)*(lvals(1)-1))*Mvals(1)+Dw(Sall2+NN(5)*(lvals(2)-1))*Mvals(2);          

    % term sum_{p=1}^3 n_p w^p_{k-n}
    DPsi4=w(Sall1)*nx(jjj)+w(Sall1+NN(4))*ny(jjj)+w(Sall1+2*NN(4))*nz(jjj);
    % the factor M_n^{l,m} - need to match dimensions
    Mnlm_expanded = Mvals(1)*(compsymvar(jj)==lvals(1))+Mvals(2)*(compsymvar(jj)==lvals(2));
    DPsi4=DPsi4.*Mnlm_expanded;
      
    j=symindex(jjj); 
    % Map to nz2-only indices
    % Find which rows in the result matrix correspond to the jj elements
    [~, row_indices] = ismember(find(jj), nz2_row_mapping);
    % Find which column in the result matrix corresponds to this j
    [~, col_index] = ismember(j, nz2_row_mapping);
    
    % add to the j-th column with factor symfactor(jjj)=\tilde{\alpha}(jjj,j)
    if col_index > 0 && any(row_indices > 0)
      valid_rows = row_indices(row_indices > 0);
      valid_mask = row_indices > 0;
      JJ(valid_rows,col_index)=JJ(valid_rows,col_index)+symfactor(jjj)*(DPsi1(valid_mask)-DPsi2(valid_mask)+DPsi3(valid_mask)-DPsi4(valid_mask));
    end
  end
  % include the factor i
  JJ=1i*JJ;
end

%%%%%%%%%%%%%%%%%%%%%% Some auxiliary functions %%%%%%%%%%%%%%%%%%

function yzext=tensorproductext(y,Dz)
% the tensorproduct used in the convection term
    szext=size(y);
    szext(1:4)=2*(szext(1:4)-1)+1;
    yzext=altzeros(szext(1:5),y(1));
    for j=1:3
        for k=1:3
              [~,convext]=convtensor(y(:,:,:,:,k),Dz(:,:,:,:,j,k));
              yzext(:,:,:,:,j)=yzext(:,:,:,:,j)+convext;
        end
    end
end

function product=P(a,b)
% elementwise product of tensors for compatibility with intlab
    product=reshape(a(:).*b(:),size(a));    
end