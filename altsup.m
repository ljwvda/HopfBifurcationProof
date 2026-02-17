function y=altsup(x)
% takes sup of interval input
% does nothing for double input

if exist('intval','file') && isintval(x(1))
    y=sup(x);
else
    y=x;
end
