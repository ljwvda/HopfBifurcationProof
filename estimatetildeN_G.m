function tildeNestimate=estimatetildeN_G(Z1target,wsol,tensors,solshapeOld,symmetry,nu,eta,Ntilde,setup)
% estimate Ntilde needed to make the Z1 tail bound smaller than Z1target

Dw=tensors.Dw;
Mw=tensors.Mw;
DMw=tensors.DMw;

% N including Z component
N=sizeshape(solshapeOld);

% New solshape
solshapeNew = solshapeOld;
solshapeNew.Nrec = [N(1) N(2) 0 0];
[nx,ny,nz,nt] = countersymtensor(solshapeNew.Nrec,solshapeNew);

% nonsymmetrized weights
weights=eta.^(abs(nx)+abs(ny)+abs(nz)+abs(nt));

% norms of (Mw)^p
Mwxi=abs(Mw).*weights;
MWxi1=Mwxi(:,:,:,:,1);
MWxi2=Mwxi(:,:,:,:,2);
MWxi3=Mwxi(:,:,:,:,3);
pMw(1)=sum(MWxi1(:));
pMw(2)=sum(MWxi2(:));
pMw(3)=sum(MWxi3(:));

% norms of D_m(Mw)
D1Mw=DMw(:,:,:,:,:,1);
D2Mw=DMw(:,:,:,:,:,2);
D3Mw=DMw(:,:,:,:,:,3);
mDMw(1)=sum(abs(D1Mw(:)).*weights(:));
mDMw(2)=sum(abs(D2Mw(:)).*weights(:));
mDMw(3)=sum(abs(D3Mw(:)).*weights(:));

% norms of D_m(w)
D1w=Dw(:,:,:,:,:,1);
D2w=Dw(:,:,:,:,:,2);
D3w=Dw(:,:,:,:,:,3);
mDw(1)=sum(abs(D1w(:)).*weights(:));
mDw(2)=sum(abs(D2w(:)).*weights(:));
mDw(3)=sum(abs(D3w(:)).*weights(:));

% norms of w^p
w=symmetrytofulltensor_ext(wsol,solshapeNew,symmetry);
wxi=abs(w).*weights;
wxi1=wxi(:,:,:,:,1);
wxi2=wxi(:,:,:,:,2);
wxi3=wxi(:,:,:,:,3);
pw(1)=sum(wxi1(:));
pw(2)=sum(wxi2(:));
pw(3)=sum(wxi3(:));

% root(Ntilde)
if exist('intval','file') && isintval(nu)
    rootN=sqrt(intval(Ntilde));
else
    rootN=sqrt(Ntilde);
end

% rescaling by appropriate factor in the bound
if exist('setup','var') && strcmp(setup,'2D')
      disp('using the 2D bound for the Z1 tail estimate')
      pMw=pMw/sqrt(nu/2);
else
      pMw=pMw/sqrt(nu/3);
end

% the three components of the tail bound
Z1tailm=max(pMw)/rootN+(3*sum(pw)/2-pw/2+mDMw+sum(mDw)-mDw)/Ntilde;
% and their maximum
% Z1tail=max(Z1tailm);     
Z1tail1=4*max(Z1tailm); 

% extra output to be able to estimate Ntilde
Z1tailcoeff=4*[max(pMw);max(3*sum(pw)/2-pw/2+mDMw+sum(mDw)-mDw)];


% solve the quadratic equation for the target value for Ntilde
tmp(2)=Z1tailcoeff(1);
tmp(1)=Z1tailcoeff(2);
tmp(3)=-Z1target;
if exist('intval','file') && isintval(nu)
    % in case this was accidentally attempted with intervals
    tmp=mid(tmp);
end
tildeNestimate=ceil(((2*tmp(1))/(-tmp(2) + sqrt(tmp(2)^2-4*tmp(1)*tmp(3))))^2);

end

