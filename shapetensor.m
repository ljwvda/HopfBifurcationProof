function ninshape = shapetensor(nx,ny,nz,nt,shape)
% codes different shapes of index sets
% SHAPE contains the information needed to define the shape

shapeType = shape.type;
if iscell(shapeType)
    shapeType = shapeType{1};
end
if isstring(shapeType)
    shapeType = char(shapeType);
end

switch shapeType
    case 'rec' % rectangle
        N=shape.Nrec;
        ninshape=(abs(nx)<=N(1) & abs(ny)<=N(2) & abs(nz)<=N(3) & abs(nt)<=N(4));
    case {'ell','ell2D'} % ellipse for Edagger and Etilde
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        ninshape=((nu*(nx.^2+ny.^2+nz.^2)).^2+(Omega*nt).^2<=tN^2);
    case 'ellHopf' % ellipse for Edagger and Etilde -- only nz=2
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        nz2 = 2*ones(size(nz));
        ninshape=((nu*(nx.^2+ny.^2+nz2.^2)).^2+(Omega*nt).^2<=tN^2);
    case 'ellHopf_nz1' % ellipse for Edagger and Etilde -- only nz=1
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        nz1 = ones(size(nz));
        ninshape=((nu*(nx.^2+ny.^2+nz1.^2)).^2+(Omega*nt).^2<=tN^2);
    case 'otherplusrec' % sum of the sets: rectangle + other shape
        N=shape.Nrec;
        anx=max(abs(nx)-N(1),0);
        any=max(abs(ny)-N(2),0);
        anz=max(abs(nz)-N(3),0);
        anz2 = 2*ones(size(anz));
        ant=max(abs(nt)-N(4),0);
        otherType = shape.other.type;
        if iscell(otherType)
            otherType = otherType{1};
        end
        if isstring(otherType)
            otherType = char(otherType);
        end
        switch otherType
            case {'ell','ell2D'}
                tN=shape.other.Nell;
                nu=shape.other.nu;
                Omega=shape.other.omega;
                ninshape=((nu*(anx.^2+any.^2+anz.^2)).^2+(Omega*ant).^2<=tN^2);
            case 'ellHopf' % Not used for making Jacobian
                tN=shape.other.Nell;
                nu=shape.other.nu;
                Omega=shape.other.omega;
                ninshape=((nu*(anx.^2+any.^2+anz2.^2)).^2+(Omega*ant).^2<=tN^2);
            case 'ellHopf_nz1' % Not used for making Jacobian
                tN=shape.other.Nell;
                nu=shape.other.nu;
                Omega=shape.other.omega;
                anz1 = ones(size(anz));
                ninshape=((nu*(anx.^2+any.^2+anz1.^2)).^2+(Omega*ant).^2<=tN^2);
            case 'rec'
                M=shape.other.Nrec;
                ninshape=(abs(anx)<=M(1) & abs(any)<=M(2) & abs(anz)<=M(3) & abs(ant)<=M(4));
            otherwise
                error('shapetensor:UnknownOtherShapeType', 'Unsupported other shape type: %s', char(otherType));
        end

    otherwise
        error('shapetensor:UnknownShapeType', 'Unsupported shape type: %s', char(shapeType));

end

end


