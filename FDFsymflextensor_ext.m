function [F,J,DFnu,Jquad,Fext,Fphase,tensors] = ...
             FDFsymflextensor_ext(x,nu,forcing,ref,solshape,symmetry,Jshape)
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

%%%%%%%%%%%%%%%%%% INDEXING %%%%%%%%%%%%%%%

if ~exist('Jshape','var')
    Jshape=solshape;
end
N=sizeshape(solshape);
M2N=max(sizeshape(Jshape),2*N);

% biggest index set needed
shapelarge.type='rec';
shapelarge.Nrec=M2N;
[nx,ny,nz,nt,ninshape,comp] = countersymtensor(M2N,shapelarge);
[symvar,symindex] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,M2N,ninshape,symmetry);
[tilden2,tilden2reci] = tildentensor(nx,ny,nz,nu);

% indices for the solution
nsol=shapetensor(nx,ny,nz,nt,solshape);
solsymvar=(symvar & nsol);
symindexsol=symindex(solsymvar);
[solnx,solny,solnz] = countersymtensor(N,solshape);
soltilden2reci=setsizetensor(tilden2reci,N);

% indices for the jacobian
njac = shapetensor(nx,ny,nz,nt,Jshape);
jacsymvar = (symvar & njac);
symindexjac = symindex(jacsymvar); 
% note that it is not necessary that nsol is a subset of njac
 
wsol=x(1:end-1);
w=symmetrytofulltensor_ext(wsol,solshape,symmetry);
w=setsizetensor(w,M2N); 
wjac=w(jacsymvar);
w=setsizetensor(w,N);

Omega=x(end);

% extend the reference for the phase condition
reflarge=altzeros([length(find(symvar)),1],nu);
reflarge(symindexsol)=ref;
refjac=reflarge(symindexjac); 

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

if nargout==7
    tensors.Dw=Dw;
    tensors.Mw=Mw;
    tensors.DMw=DMw;
end

% calculate F
forcingsol=fulltosymmetrytensor_ext(setsizetensor(forcing,N),solshape,symmetry);
F1 = (1i*Omega*nt(solsymvar)+nu*tilden2(solsymvar)).*wsol - forcingsol;
F1large2N=setsizetensor(symmetrytofulltensor_ext(F1,solshape,symmetry),2*N);
F2large2N=tensorproductext(Mw,Dw); 
F3large2N=tensorproductext(w,DMw); 
Fext=F1large2N+1i*(F2large2N-F3large2N); % this has tensor data structure

Fsol=fulltosymmetrytensor_ext(setsizetensor(Fext,N),solshape,symmetry);

% phase condition
Fphase=1i*ref'*(nt(solsymvar).*wsol);

F=[Fsol;Fphase];

if nargout==1 % no Jacobian computed if not needed
    return
end


%%%%%%%%%% JACOBIAN %%%%%%%%%%%% 

% derivative of the nonlinear part in omega
Jnonlinear=Dnonlinear(w,Dw,Mw,DMw,N,Jshape,symmetry);

MJ=length(find(jacsymvar));
J=altzeros([MJ+1,MJ+1],nu);
J(1:MJ,1:MJ)=Jnonlinear; 

diagMJ=(speye(MJ+1)==1);
diagMJ(end,end)=false;

% derivative of the time derivative term w.r.t. Omega
J(1:end-1,end)=1i*nt(jacsymvar).*wjac; 

% derivative of the time derivative term w.r.t. omega
J(diagMJ)=J(diagMJ)+1i*Omega*nt(jacsymvar);
% Jquad is an extra output containing the derivative of all quadratic terms
Jquad=J; 
% derivative of the Laplacian term
J(diagMJ)=J(diagMJ)+nu*tilden2(jacsymvar);

% derivative of phase condition
J(end,1:end-1)=1i*(nt(jacsymvar).*refjac)';

% DFnu is an extra output containing the derivative w.r.t. the parameter nu
DFnu=[tilden2(jacsymvar).*wjac;0];

end  % of main function 

function JJ=Dnonlinear(w,Dw,Mw,DMw,N,Jshape,symmetry)
% determines the Jacobian of the nonlinear part

  % initialize indices
  M=sizeshape_Hopf(Jshape); 
  MN=max(M,N);
  [nx,ny,nz,nt,ninshape,comp] = countersymtensor(MN,Jshape);
  [symvar,symindex,symfactor] = symmetryindicestensor_ext(nx,ny,nz,nt,comp,MN,ninshape,symmetry);
  [~,tilden2reci] = tildentensor(nx,ny,nz,Dw(1));
  JM=length(find(symvar));
  JJ=altzeros([JM,JM],Dw(1));
  [nsymvar,NN]=makeindices(N,nx,ny,nz,nt,symvar);
  compsymvar=comp(symvar);
  compind{1}=(compsymvar==1);
  compind{2}=(compsymvar==2);
  compind{3}=(compsymvar==3);

  % for loop to fill the columns of the matrix 
  for jjj=find(symindex)'      
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
    % will lead to a nonvanishing contribution
    jj=(abs(nsx)<=N(1) & abs(nsy)<=N(2) & abs(nsz)<=N(3) & abs(nst)<=N(4));
    
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
    % the factor M_n^{l,m}
    DPsi4=DPsi4.*Mnlm(jj);
      
    j=symindex(jjj); 
    % add to the j-th column with factor symfactor(j)=\tilde{\alpha}(jj,j)
    if any(jj)
      JJ(jj,j)=JJ(jj,j)+symfactor(jjj)*(DPsi1-DPsi2+DPsi3-DPsi4);
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

