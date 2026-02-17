function [tilden2,tilden2reci] = tildentensor(nx,ny,nz,intvalcheck)
% contructs \tilde{n}^2 and its reciprocal

tilden2=nx.^2+ny.^2+nz.^2;
if exist('intval','file') && isintval(intvalcheck)
  tilden2reci=reshape(1./intval(tilden2(:)),size(tilden2));
  tilden2reci(tilden2==0)=0;
else
  tilden2reci=1./tilden2;  
  tilden2reci(tilden2==0)=0;
end

end

