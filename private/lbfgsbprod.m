function p = lbfgsbprod(H,g)
%lbfgshprod  computes products with the L-BFGS matrix B.
%
%   p = lbfgsbprod(H,g)  returns  p = B*g.
%
%   The product is computed using the factoriation
%   [(9.15),p.231] described in Nocedal and Wright, 1999.

%   See also lbfgsadd, lbfgsdel, lbfgsupdate, lbfgsinit.

% $Id$

% ----------------------------------------------------------------------
% Explicit matrix formulation
% ----------------------------------------------------------------------

%jMax = H.jMax;
%jNew = H.jNew;
%jOld = H.jOld;

% |---|---|---|---|---|
%       ^   ^
%       |   |
%      old new
%   4   5   1   2   3    1=newest, 2=2nd newest,... 5=oldest

%if false
%   Sk = H.S(:,[jNew-1:-1:1,jMax:-1:jNew]);
%   Yk = H.Y(:,[jNew-1:-1:1,jMax:-1:jNew]);
%   Lk = (Sk' * Yk) .* (repmat((1:jMax)',1,jMax) > repmat((1:jMax),jMax,1)) ;
%   Dk = diag(sum(Sk.*Yk));

%   deltak = H.delta;

%   M = [deltak*Sk'*Sk, Lk; Lk', -Dk];

%   p = [deltak*Sk'; Yk'] * g;
%   p = M \ p;
%   p = [deltak*Sk, Yk] * p;
%   p = deltak * g - p;
%end

% ----------------------------------------------------------------------
% Formulation directly based on arrays in permuted form :-)
% ----------------------------------------------------------------------

if ~isempty(H.ML)

   % This code works much faster for larger problems
   v1 = H.delta * (g' * H.S)';
   v2 = (g' * H.Y)';
   p  = [v1(H.valid); v2(H.valid)];

%  v1 = H.delta * (g' * H.S(:,H.valid))';
%  v2 = (g' * H.Y(:,H.valid))';
%  q = [v1; v2];

   p = H.MU \ (H.ML \ p); % H.M \ p (optionally: use linsolve)

%  v1 = H.S(:,H.valid) * p(1:H.rank) * H.delta;
%  v2 = H.Y(:,H.valid) * p(H.rank+1:2*H.rank); 
%  q = v1 + v2;

   % This code works much faster for larger problems
   pe = zeros(size(H.S,2),1); pe(H.valid) = p(1:H.rank);
   v1 = H.S * (pe * H.delta);
   pe = zeros(size(H.Y,2),1); pe(H.valid) = p(H.rank+1:2*H.rank);
   v2 = H.Y * pe; 
   p  = v1 + v2;
else
   p = 0;
end

p = H.delta * g - p;