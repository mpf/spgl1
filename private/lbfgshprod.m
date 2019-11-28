function p = lbfgshprod(H,g)
%lbfgshprod  computes products with the inverse of the L-BFGS matrix H.
%
%   p = lbfgshprod(H,g)  returns  p = H*g.
%
%   The product is computed using the recursive algorithm (Algorithm
%   9.2, p.226) described in Nocedal and Wright, 1999.

%   See also lbfgsadd, lbfgsdel, lbfgsupdate, lbfgsinit.

% $Id: lbfgshprod.m 385 2008-09-02 22:22:41Z mpf $

jMax = H.jMax;
jNew = H.jNew;
alfa = zeros(jMax,1);

% |---|---|---|---|---|
%       ^   ^
%       |   |
%      old new
%   4   5   1   2   3    1=newest, 2=2nd newest,... 5=oldest

p = g;
for i = 1:jMax                % Run from "newest" to "oldest"
    j       = jNew + i - 1;
    k       = mod(j - 1,jMax) + 1;
    s       = H.S(:,k);
    y       = H.Y(:,k);
    rho     = H.r(  k);
    alfa(k) = rho*(s'*p);
    p       = p - alfa(k)*y;
end
p = H.gamma * p; % H0 = H.gamma * Id
for i = jMax:-1:1             % Run from "oldest" to "newest"
    j       = jNew + i - 1;
    k       = mod(j - 1,jMax) + 1;
    s       = H.S(:,k);
    y       = H.Y(:,k);
    rho     = H.r(  k);
    beta    = rho*(y'*p);
    p       = p + (alfa(k) - beta)*s;
end



