function wout = setsizetensor(win,N)
% sets the size of the (full) tensor to N
% pads with zeros in dimensions that are too small
% truncates in dimensions that are too large 

sz=size(win);
M=(sz(1:4)-1)/2;

otherdim=prod(sz(5:end));
win=reshape(win,[sz(1:4),otherdim]);

L=max(N,M);
wlarge=altzeros([2*L+1,otherdim],win(1));
wlarge(L(1)+1+(-M(1):M(1)),L(2)+1+(-M(2):M(2)),...
                     L(3)+1+(-M(3):M(3)),L(4)+1+(-M(4):M(4)),:)=win;
wout=wlarge(L(1)+1+(-N(1):N(1)),L(2)+1+(-N(2):N(2)),...
                     L(3)+1+(-N(3):N(3)),L(4)+1+(-N(4):N(4)),:);
wout=reshape(wout,[2*N+1,sz(5:end)]);                 

end

