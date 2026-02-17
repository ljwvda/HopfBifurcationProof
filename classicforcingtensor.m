function forcing = classicforcingtensor
% returns the classical value of the forcing f^omega
% as a small tensor

forcing=zeros([3,3,1,1,3]);
forcing(3,3,1,1,3)=-1;
forcing(1,3,1,1,3)=1;
forcing(3,1,1,1,3)=1;
forcing(1,1,1,1,3)=-1;

end

