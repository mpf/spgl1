function H = lbfgsinit(n,k,dscale)
%lbfgsinit  initializes an L-BFGS matrix.
%
%   H = lbfgsinit(n,k,dscale)
%
%   See also lbfgsadd, lbfgsdel, lbfgshprod, lbfgsupdate.

% $Id: lbfgsinit.m 385 2008-09-02 22:22:41Z mpf $

if nargin < 3 || isempty(dscale)
   dscale = 1;
end

H.jNew  = 1;
H.jOld  = 1;
H.jMax  = k;
H.S     = zeros(n,k);
H.Y     = zeros(n,k);

% Data specific to H
H.gamma = dscale;
H.r     = zeros(1,k);

% Data specific to B
H.delta = 1 / dscale;
H.valid = false(k,1);
H.rank  = 0;
H.STS   = zeros(k,k);
H.L     = zeros(k,k);
H.D     = zeros(k,k);
H.M     = zeros(2*k,2*k);
H.ML    = [];
H.MU    = [];

% Status information and options
H.status = 0; % 0=updated, 1=corrected, 2=no update