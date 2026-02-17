function N=sizeshape(shape)
% determines the tightest tensor that contains the shape

switch shape.type
    case {'rec'} % rectangle
        N=shape.Nrec;
    case {'ell'} % for Edagger (for Adagger) and Etilde
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        N(1:3)=ceil(altsup(sqrt(tN/nu)));
        N(4)=ceil(altsup(tN/Omega));
    case {'ell2D'} % ellipse for 2D
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        N(1:2)=ceil(altsup(sqrt(tN/nu)));
        N(3)=0; % this saves quite a bit of memory
        N(4)=ceil(altsup(tN/Omega));
    case {'ellHopf'} % ellipse for 2D with nz=2 only
        tN=shape.Nell;
        nu=shape.nu;
        Omega=shape.omega;
        N(1:2)=ceil(altsup(sqrt(tN/nu)));
        N(3)=0; % this saves quite a bit of memory
        % N(4)=ceil(altsup(tN/Omega));
        N(4)=0; % Adjusted
    case {'otherplusrec'} % sum of the sets: rectangle + other shape
        N=shape.Nrec+sizeshape(shape.other);
end

