function [wout,wfull] = conjugatesymmetric_ext(win,shape,symmetry)
% returns conjugate symmetric approximation of the input
% works for both intvals and floats

wfull = symmetrytofulltensor_ext(win,shape,symmetry);
wfull = (wfull+conj(wfull(end:-1:1,end:-1:1,end:-1:1,end:-1:1,:)))/2;
wout = fulltosymmetrytensor_ext(wfull,shape,symmetry);

end

