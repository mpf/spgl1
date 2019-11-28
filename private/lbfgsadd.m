function H = lbfgsadd(H,k,dscale)
%lbfgsadd  adds k variables to the L-BFGS matrix H.
%
%   H = lbfgsadd(H,k,dscale) expands H to accomodate k new variables.
%   The new diagonal part of H (H0) gets set to dscale.
%
%   See also lbfgsdel, lbfgshprod, lbfgsupdate, lbfgsinit.

% $Id: lbfgsadd.m 385 2008-09-02 22:22:41Z mpf $

jMax = H.jMax;

H.S  = [ H.S ;    zeros(k,jMax)   ];
H.Y  = [ H.Y ;    zeros(k,jMax)   ];
H.H0 = [ H.H0; repmat(dscale,k,1) ];
