function [H,noup] = lbfgsupdate(H, step, p, g1, g2)
%lbfgsupdate  updates the L-BFGS matrix.
%
%   [H,noup] = lbfgsupdate(H,step,p,g1,g2) updates H.  The flag noup
%   indicates if the update was skipped.
%
%   See also lbfgsadd, lbfgsdel, lbfgshprod, lbfgsinit.

% 02 Sept 2008: Continue to "age" the vectors even if an update is skipped.
%               This means that the pointers will get advanced and that the
%               oldest vectors won't participate in the next Hprod.

% Storage scheme: jOld points to the location with the oldest
% vectors. The next location has the newest vectors.  For examples, if
% jMax = 5, then
%
% |---|---|---|---|---|
%       ^   ^
%       |   |
%      old new
%   4   5   1   2   3    1=newest, 2=2nd newest,...,5=oldest
%
% $Id: lbfgsupdate.m 389 2008-09-03 06:03:53Z mpf $

% Get size of storage and pointer to the last updated vectors.
jMax = H.jMax;
jOld = H.jOld;

gtp1 = g1'*p;
gtp2 = g2'*p;
noup = gtp2 <= 0.91*gtp1; % Curvature requirement

if ~noup
   H.status = 0; % Update
  
   % See Nocedal & Wright, Second Edition, page 537
   s   = step * p;
   y   = g2 - g1;
   bs  = lbfgsbprod(H,s);
   sbs = s' * bs;
   yts = step * (gtp2 - gtp1); % TODO: y'*s?

   if yts < 0.2 * sbs
      theta = (0.8 * sbs) / (sbs - yts);
      y     = theta * y   + (1-theta)*bs;
      yts   = theta * yts + (1-theta)*sbs;
      H.status = 1; % Correction
   end
end


if noup
   % Shift all vectors and update factorization;
   % Do not add any new vectors.

   H.valid(jOld) = false;
   H.r(jOld) = 0;    % Causes the corresponding vectors to not participate
                     % when computing subsequent products with H. 
   H.status = 2; % No update;
else
%   y   = g2 - g1;
%   yts = step*(gtp2 - gtp1);
   yty = y'*y;

   % Replace oldest vectors in S and Y.
   H.S(:,jOld) = s; % TODO: step*p;
   H.Y(:,jOld) = y;

   % Update initial matrix H0 = H.gamma * Id
   H.gamma   = yts / yty; % s'y / y'y
   H.r(jOld) = 1 / yts;
   
   % Update data for B
   sTS = step*p'*H.S;
   H.delta = yty / yts; % y'y / s'y
   H.valid(jOld)  = true;
   H.STS(jOld,:)  = sTS;
   H.STS(:,jOld)  = sTS';
   H.D(jOld,jOld) = yts;
   H.L(jOld,:)    = s'*H.Y; % TODO: was step*p'
   H.L(:,jOld)    = 0;
end

% Update and factorize matrix M for B
H.rank      = sum(H.valid);
H.M         = [H.delta * H.STS, H.L; H.L', -H.D];
[H.ML,H.MU] = lu(H.M([H.valid;H.valid],[H.valid;H.valid]));

% Advance pointers: the newest vectors are now where the oldest used to be.
% The oldest vectors are always one spot "behind" the newest.
H.jNew = jOld;
if jOld == 1
   H.jOld = jMax;
else
   H.jOld = jOld - 1;
end



% ---------------------
% DEBUG - DEBUG - DEBUG
% ---------------------

%disp(H.rank)
%g0 = randn(size(H.S,1),1);
%g1 = lbfgsbprod(H,g0);
%g2 = lbfgshprod(H,g1);
%disp(norm(g2-g0))
%plot(g0); hold on; plot(g2,'.'); hold off; pause;
