function H = lbfgsdel(H,kdel)
%lbfgsdel  deletes variables from an L-BFGS matrix.
%
%   H = lbfgsdel(H,kdel) wipes out all information in H about the
%   variables indexed by kdel.
%
%   See also lbfgsadd, lbfgshprod, lbfgsupdate, lbfgsinit.

% $Id: lbfgsdel.m 385 2008-09-02 22:22:41Z mpf $

jMax = H.jMax;

H.S (kdel,:) = [];
H.Y (kdel,:) = [];
H.H0(kdel  ) = [];

% The inner products s'y may have changed. Recompute these and store
% their inverses in r.  If a pair (s,y) is nearly orthogonal, set r=0;
% this causes them not to participate in products with H.
for j = 1:jMax
    s   = H.S(:,j);
    y   = H.Y(:,j);
    sty = s'*y;
    if sty > 1e-5*norm(s)*norm(y)
       H.r(j) = 1/sty;
    else
       H.r(j) = 0;
    end
end