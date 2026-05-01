function N=sizeshape_Hopf(shape)
% determines the tightest tensor that contains the shape
% Adjusted version of sizeshape.m for Hopf bifurcation

% N = zeros(1,4);
% 
% shapeType = shape.type;
% if iscell(shapeType)
%     shapeType = shapeType{1};
% end
% if isstring(shapeType)
%     shapeType = char(shapeType);
% end

switch shape.type
    case 'rec' % rectangle
        N=shape.Nrec;
    case 'ell' % for Edagger (for Adagger) and Etilde
        tN=shape.Nell;
        nu=shape.nu;
        N(1:3)=ceil(altsup(sqrt(tN/nu)));
        %N(4)=ceil(altsup(tN/Omega));
        N(4)=0; % adjusted
    case 'ell2D' % ellipse for 2D
        tN=shape.Nell;
        nu=shape.nu;
        N(1:2)=ceil(altsup(sqrt(tN/nu)));
        N(3)=0; % this saves quite a bit of memory
        % N(4)=ceil(altsup(tN/Omega));
        N(4)=0; % Adjusted
    case {'ellHopf','ellHopf_nz1'} % ellipse for 2D with fixed nz
        tN=shape.Nell;
        nu=shape.nu;
        N(1:2)=ceil(altsup(sqrt(tN/nu)));
        N(3)=0; % this saves quite a bit of memory
        % N(4)=ceil(altsup(tN/Omega));
        N(4)=0; % Adjusted
    case 'otherplusrec' % sum of the sets: rectangle + other shape
        N=shape.Nrec+sizeshape_Hopf(shape.other);
    otherwise
        error('sizeshape_Hopf:UnknownShapeType', 'Unsupported shape type: %s', char(shapeType));
end

end

