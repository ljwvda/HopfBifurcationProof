function [symvar,symindex,symfactor,grouporder,multiplicity,symindexextra] ...
    = symmetryindicestensor_ext(nx,ny,nz,nt,comp,N,ninshape,symmetry)
% creates several indices and auxiliary variables related to the symmetry
% symvar: indicates which are the symmetry reduced indices (Jred)
% symindex: for each index in Jsym tells us the (unique) element of Jred in its orbit
% symfactor: contains the values of \tilde{\alpha} for indices in Jsym
% grouporder: is the order of the group
% multiplicity: is the order of the stabilizer subgroup (for indices in Jsym)
% symindexextra: for each index in Jred tells us which elements are in its orbit

% Extended with cases for second and third Hopf bifurcation point

% determines a fundamental domain
switch symmetry
    case 0 % no symmetries
        nsym=true(size(nx));
        grouporder=1;
    case 1 % Sz symmetry
        nsym=(nz>=0);
        grouporder=2;
    case 11 % Sz symmetry - for ellHopf case
        nsym=(nz>=0);
        grouporder=2;
    case 2 % Sz and DSx and P2D symmetry
        nsym=(nz>=0 & nx>=0);
        grouporder=8;
    case 4 % 2D flow and SxSy and DSx and P4SxR and P2D symmetry
        nsym=(nx>=0 & ny >=0 & ny <= nx);
        grouporder=16;
        independentofz=(nz==0);
        zcomponly=(comp==3);
        nsym= (nsym & independentofz & zcomponly);
    case 25 % Sz and SxSy and DSx symmetry
        nsym=(nz>=0 & nx>=0 & ny>=0);
        grouporder=8;
    case 26 % Sz and SxSy and DSx symmetry - for ellHopf case
        nsym=(nz>=0 & nx>=0 & ny>=0);
        grouporder=8;
end

% "average" over the group action (without normalizing with the group order)
symfactor=symmetrizetensor_ext(nsym,symmetry);

ncenter=(nx==0 & ny==0 & nz==0 & nt==0);
% symvar=(symfactor~=0 & nsym & ninshape & ~ncenter);
% For nz=2 case, we do take the center (0,0,2) into account
if (symmetry == 11 || symmetry==26)
    symvar=(symfactor~=0 & nsym & ninshape);
else
    symvar=(symfactor~=0 & nsym & ninshape & ~ncenter);
end

% order of the stabilizer for symmetry indices
multiplicity=abs(symfactor);
multiplicity(multiplicity==0)=1; % prevent 0/0 later;
% symfactor is \tilde{\alpha} for symmetry indices
symfactor=symfactor./multiplicity;

% ordering of the symmetry reduced indices
numbersymvar=length(find(symvar));
symvarindex=(1:numbersymvar)';
symredindex=zeros([2*N+1,3]);
symredindex(symvar)=symvarindex;

% use the group action to determine the orbits of the symmetry reduced indices
wtmp{1}=symredindex;
% tracks the group action (beta) on the indices
switch symmetry
    case 0 % no symmetry
        grouporder=1;
    case 1 % Sz symmetry
        wtmp{2}=symredindex(:,:,end:-1:1,:,:);   %Sz symmetry
        grouporder=2;
    case 11 % Sz symmetry for ellHopf case
        wtmp{2}=symredindex(:,:,end:-1:1,:,:);   %Sz symmetry
        grouporder=2;
    case 2 % Sz and DSx and P2D symmetry
        wtmp{2}=symredindex(:,:,end:-1:1,:,:);   %Sz symmetry
        for k=3:4
            wtmp{k}=wtmp{k-2}(end:-1:1,:,:,:,:);   %DSx symmetry
        end
        for k=5:8
            wtmp{k}=wtmp{k-4};  %P2D symmetry
        end
        grouporder=8;
    case 4  % 2D flow
        wtmp{2}=symredindex(end:-1:1,end:-1:1,:,:,:);   %SxSy symmetry
        for k=3:4
            wtmp{k}=wtmp{k-2}(end:-1:1,:,:,:,:);   %DSx symmetry
        end
        for k=5:8
            wtmp{k}=permute(wtmp{k-4},[2 1 3 4 5]);  %P4SxR symmetry
            wtmp{k}(:,:,:,:,[1 2])=wtmp{k}(:,:,:,:,[2 1]);
        end
        for k=9:16
            wtmp{k}=wtmp{k-8};  %P2D symmetry
        end
        grouporder=8;
    case 25 % Sz and SxSy and DSx symmetry 
        wtmp{2}=symredindex(:,:,end:-1:1,:,:);   %Sz symmetry
        for k=3:4
            wtmp{k}=wtmp{k-2}(end:-1:1,end:-1:1,:,:,:);   %SxSy symmetry
        end
        for k=5:8
            wtmp{k}=wtmp{k-4}(end:-1:1,:,:,:,:);   %DSx symmetry
        end
        grouporder=8;
    case 26 % Sz and SxSy and DSx symmetry - for ellHopf case
        wtmp{2}=symredindex(:,:,end:-1:1,:,:);   %Sz symmetry
        for k=3:4
            wtmp{k}=wtmp{k-2}(end:-1:1,end:-1:1,:,:,:);   %SxSy symmetry
        end
        for k=5:8
            wtmp{k}=wtmp{k-4}(end:-1:1,:,:,:,:);   %DSx symmetry
        end
        grouporder=8;
end

% symindex: for each index in Jsym tells us the (unique) element of Jred in its orbit
symindex=symredindex;
for k=2:grouporder
    nonzero=wtmp{k}>0;
    symindex(nonzero)=wtmp{k}(nonzero);
end

if nargout==6
    % symindexextra: for each index in Jred tells us which elements are in its orbit
    symindexextra=zeros(nnz(symredindex),grouporder);
    for k=1:grouporder
        symindexextra(:,k)=sorted(wtmp{k});
    end
end

end

function wsort=sorted(win)
% finds the indices in J where the indices in Jred got mapped to
% sorted according to the ordering of Jred
pos=find(win>0);
val=win(pos);
[~,ind]=sort(val);
wsort=pos(ind);
end
