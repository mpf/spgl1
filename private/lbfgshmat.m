function M = lbfgshmat(H)
%lbfgshmat  forms explicit inverse of the L-BFGS matrix, H.
%
%   M = lbfgshprod(H)  returns  L-BFGS inverse matrix.
%
% This routine is for debugging purposes only.

%   See also lbfgshprod, lbfgsinit.

% $Id$

jMax = H.jMax;
jNew = H.jNew;
alfa = zeros(jMax,size(H.S,1));

% |---|---|---|---|---|
%       ^   ^
%       |   |
%      old new
%   4   5   1   2   3    1=newest, 2=2nd newest,... 5=oldest

M = eye(size(H.S,1));

for i = 1:jMax                % Run from "newest" to "oldest"
    j       = jNew + i - 1;
    k       = mod(j - 1,jMax) + 1;
    s       = H.S(:,k);
    y       = H.Y(:,k);
    rho     = H.r(  k);
    alfa(k,:) = rho*(s'*M);
    M       = M - y*alfa(k,:);
end
M = H.gamma * M; % H0 = H.gamma * Id
for i = jMax:-1:1             % Run from "oldest" to "newest"
    j       = jNew + i - 1;
    k       = mod(j - 1,jMax) + 1;
    s       = H.S(:,k);
    y       = H.Y(:,k);
    rho     = H.r(  k);
    beta    = rho*(y'*M);
    M       = M + s*(alfa(k,:) - beta);
end



