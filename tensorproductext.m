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