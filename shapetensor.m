function ninshape = shapetensor(nx,ny,nz,nt,shape)
% codes different shapes of index sets
% SHAPE contains the information needed to define the shape

switch shape.type
    case {'rec'} % rectangle
        N=shape.Nrec;
        ninshape=(abs(nx)<=N(1) & abs(ny)<=N(2) & abs(nz)<=N(3) & abs(nt)<=N(4));
    case {'ell','ell2D'} % ellipse for Edagger and Etilde
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        ninshape=((nu*(nx.^2+ny.^2+nz.^2)).^2+(Omega*nt).^2<=tN^2);
    case {'ellHopf'} % ellipse for Edagger and Etilde -- only nz=2
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        nz2 = 2*ones(size(nz));
        ninshape=((nu*(nx.^2+ny.^2+nz2.^2)).^2+(Omega*nt).^2<=tN^2);
    case 'otherplusrec' % sum of the sets: rectangle + other shape
        N=shape.Nrec;
        anx=max(abs(nx)-N(1),0);
        any=max(abs(ny)-N(2),0);
        anz=max(abs(nz)-N(3),0);
        anz2 = 2*ones(size(anz));
        ant=max(abs(nt)-N(4),0);
        switch shape.other.type
            case {'ell','ell2D'}
                tN=shape.other.Nell;
                nu=shape.other.nu;
                Omega=shape.other.omega;
                ninshape=((nu*(anx.^2+any.^2+anz.^2)).^2+(Omega*ant).^2<=tN^2);
            case {'ellHopf'} % Not used for making Jacobian
                tN=shape.other.Nell;
                nu=shape.other.nu;
                Omega=shape.other.omega;
                ninshape=((nu*(anx.^2+any.^2+anz2.^2)).^2+(Omega*ant).^2<=tN^2);
            case {'rec'}
                M=shape.other.Nrec;
                ninshape=(abs(anx)<=M(1) & abs(any)<=M(2) & abs(anz)<=M(3) & abs(ant)<=M(4));
        end

end


