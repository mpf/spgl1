function [x,r,g,info] = spgl1( A, b, tau, sigma, x, options )
%SPGL1  Solve basis pursuit, basis pursuit denoise, and LASSO
%
% [x, r, g, info] = spgl1(A, b, tau, sigma, x0, options)
%
% ---------------------------------------------------------------------
% Solve the basis pursuit denoise (BPDN) problem
%
% (BPDN)   minimize  ||x||_1  subj to  ||Ax-b||_2 <= sigma,
%
% or the l1-regularized least-squares problem
%
% (LASSO)  minimize  ||Ax-b||_2  subj to  ||x||_1 <= tau.
%
% Setting the optional parameter mu uses updated A and b:
% A = [A; sqrt(\mu)I] and b = [b; 0].
% ---------------------------------------------------------------------
%
% INPUTS
% ======
% A        is an m-by-n matrix, explicit or an operator.
%          If A is a function, then it must have the signature
%
%          y = A(x,mode)   if mode == 1 then y = A x  (y is m-by-1);
%                          if mode == 2 then y = A'x  (y is n-by-1).
%
% b        is an m-vector.
% tau      is a nonnegative scalar; see (LASSO).
% sigma    if sigma != inf or != [], then spgl1 will launch into a
%          root-finding mode to find the tau above that solves (BPDN).
%          In this case, it's STRONGLY recommended that tau = 0.
% x0       is an n-vector estimate of the solution (possibly all
%          zeros). If x0 = [], then SPGL1 determines the length n via
%          n = length( A'b ) and sets  x0 = zeros(n,1).
% options  is a structure of options from spgSetParms. Any unset options
%          are set to their default value; set options=[] to use all
%          default values.
%
% OUTPUTS
% =======
% x        is a solution of the problem
% r        is the residual, r = b - Ax
% g        is the gradient, g = -A'r
% info     is a structure with the following information:
%          .tau     final value of tau (see sigma above)
%          .rNorm   two-norm of the optimal residual
%          .rGap    relative duality gap (an optimality measure)
%          .gNorm   Lagrange multiplier of (LASSO)
%          .stat    = 1 found a BPDN solution
%                   = 2 found a BP sol'n; exit based on small gradient
%                   = 3 found a BP sol'n; exit based on small residual
%                   = 4 found a LASSO solution
%                   = 5 error: too many iterations
%                   = 6 error: linesearch failed
%                   = 7 error: found suboptimal BP solution
%                   = 8 error: too many matrix-vector products
%          .time    total solution time (seconds)
%          .nProdA  number of multiplications with A
%          .nProdAt number of multiplications with A'
%
% OPTIONS
% =======
% Use the options structure to control various aspects of the algorithm:
%
% options.fid         File ID to direct log output
%        .verbosity   0=quiet, 1=some output, 2=more output.
%        .iterations  Max. number of iterations (default if 10*m).
%        .bpTol       Tolerance for identifying a basis pursuit solution.
%        .optTol      Optimality tolerance (default is 1e-4).
%        .decTol      Larger decTol means more frequent Newton updates.
%        .subspaceMin 0=no subspace minimization, 1=subspace minimization.
%
% EXAMPLE
% =======
%   m = 120; n = 512; k = 20; % m rows, n cols, k nonzeros.
%   p = randperm(n); x0 = zeros(n,1); x0(p(1:k)) = sign(randn(k,1));
%   A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
%   b  = A*x0 + 0.005 * randn(m,1);
%   opts = spgSetParms('optTol',1e-4);
%   [x,r,g,info] = spgl1(A, b, 0, 1e-3, [], opts); % Find BP sol'n.
%
% AUTHORS
% =======
%  Ewout van den Berg (vandenberg.ewout@gmail.com)
%  Michael P. Friedlander (mpf@math.ucdavis.edu)
%
% BUGS
% ====
% Please send bug reports or comments to
%            Michael P. Friedlander (mpf@math.ucdavis.edu)
%            Ewout van den Berg (vandenberg.ewout@gmail.com)

% 15 Apr 07: First version derived from spg.m.
% 17 Apr 07: Added root-finding code.
% 18 Apr 07: sigma was being compared to 1/2 r'r, rather than
%            norm(r), as advertised.  Now immediately change sigma to
%            (1/2)sigma^2, and changed log output accordingly.
% 24 Apr 07: Added quadratic root-finding code as an option.
% 24 Apr 07: Exit conditions need to guard against small ||r||
%            (ie, a BP solution).  Added test1,test2,test3 below.
% 15 May 07: Trigger to update tau is now based on relative difference
%            in objective between consecutive iterations.
% 15 Jul 07: Added code to allow a limited number of line-search
%            errors. 
% 23 Feb 08: Fixed bug in one-norm projection using weights. Thanks
%            to Xiangrui Meng for reporting this bug.
% 26 May 08: The simple call spgl1(A,b) now solves (BPDN) with sigma=0.
% 18 Mar 13: Reset f = fOld if curvilinear line-search fails.
%            Avoid computing the Barzilai-Borwein scaling parameter
%            when both line-search algorithms failed.
% 09 Sep 13: Recompute function information at new x when tau decreases.
%            Fixed bug in subspace minimization. Thanks to Yang Lei
%            for reporting this bug.

% Version 2
% - history flag (false by default)
% - c, mu
% - Keep track of the best dual value
% - Keep track of the best tau'

%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected-Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------
REVISION = '2.0';
DATE     = '29 Apr 2015';

tic; % Start your watches!

%----------------------------------------------------------------------
% Check arguments and determine solver mode
%----------------------------------------------------------------------
if ~exist('options','var'), options = []; end
if ~exist('x','var'), x = []; end
if ~exist('sigma','var'), sigma = []; end
if ~exist('tau','var'), tau = []; end

if nargin < 2 || isempty(b) || isempty(A)
   error('At least two arguments are required');
elseif (isempty(sigma) && ~isempty(tau))   
   singleTau = true;
else
   if (isempty(tau)), tau = 0; end;
   if (isempty(sigma)), sigma = 0; end;
   singleTau = false;
end


%----------------------------------------------------------------------
% Get problem size and parameters
%----------------------------------------------------------------------
m = length(b);

% Determine initial x, vector length n, and see if problem is complex
explicit = ~(isa(A,'function_handle'));
if isnumeric(A)
   n = size(A,2);
   realx = isreal(A) && isreal(b);
else
   if isempty(x)
      x = Aprod(b,2);
   end
   n = length(x);
   realx = isreal(x) && isreal(b);
end
if isempty(x)
   x = zeros(n,1);
end


%----------------------------------------------------------------------
% Apply default parameters where needed
%----------------------------------------------------------------------
defaultopts = spgSetParms();
defaultopts.iterations = 10*m;
options = spgSetParms(defaultopts, options);

if (~isnan(options.iscomplex)), realx = (options.iscomplex == 0); end


%----------------------------------------------------------------------
% Check if all weights (if any) are strictly positive. In previous
% versions we also checked if the number of weights was equal to n.
% In the case of multiple measurement vectors, this no longer needs
% to apply, so the check was removed.
%----------------------------------------------------------------------
if (~isempty(options.weights))
   if (any(~isfinite(options.weights)))
      error('Entries in options.weights must be finite');
   end
   if (any(options.weights <= 0))
      error('Entries in options.weights must be strictly positive');
   end
else
   options.weights = 1;
end


%----------------------------------------------------------------------
% Initialize the local variables
%----------------------------------------------------------------------

% Problem parameters
bNorm         = norm(b,2);
mu            = options.mu;
weights       = options.weights;

% Operational variables
stat          = 0;
nPrevVals     = options.nPrevVals;
lastFv        = -inf(nPrevVals,1);  % Last m function values.
stepG         = 1;                  % Step length for projected gradient.

% Compute modes
rfModePrimal  = (options.rootfindMode == 0);
rfTolDual     = options.rootfindTol;
hybridMode    = options.hybridMode;

% Bounds
stepMin       = options.stepMin;
stepMax       = options.stepMax;
maxIts        = options.iterations;
maxMatvec     = max(3,options.maxMatvec);
maxLineErrors = 10;     % Maximum number of line-search failures.

% Tolerance levels
bpTol         = options.bpTol;
lsTol         = options.lsTol;
optTol        = options.optTol;
decTol        = options.decTol;
relgapMinF    = options.relgapMinF; % gap / max(relgapMinF, f)
relgapMinR    = options.relgapMinR; % rNorm / max(relgapMinR, rNorm)
rootfindTol   = options.rootfindTol; % (dual - sigma) / (primal - sigmal)

% Runtime statistics
iter          = 0;
nProdA        = 0;
nProdAt       = 0;
nLineTot      = 0; % Total no. of linesearch steps.
nNewton       = 0;
timeProject   = 0;
timeMatProd   = 0;

% Output configuration
fid           = options.fid;
logLevel      = options.verbosity;
printTau      = false;
history       = options.history;


%----------------------------------------------------------------------
% Exit condition constants
%----------------------------------------------------------------------

EXIT_ROOT_FOUND    = 1;
EXIT_BPSOL_FOUND   = 2;
EXIT_LEAST_SQUARES = 3;
EXIT_OPTIMAL       = 4;
EXIT_ITERATIONS    = 5;
EXIT_LINE_ERROR    = 6;
EXIT_SUBOPTIMAL_BP = 7;
EXIT_MATVEC_LIMIT  = 8;

EXIT_MESSAGES    = cell(8,1);
EXIT_MESSAGES{1} = 'EXIT -- Found a root';
EXIT_MESSAGES{2} = 'EXIT -- Found a BP solution';
EXIT_MESSAGES{3} = 'EXIT -- Found a least-squares solution';
EXIT_MESSAGES{4} = 'EXIT -- Optimal solution found';
EXIT_MESSAGES{5} = 'ERROR EXIT -- Too many iterations';
EXIT_MESSAGES{6} = 'ERROR EXIT -- Linesearch error';
EXIT_MESSAGES{7} = 'EXIT -- Found a suboptimal BP solution';
EXIT_MESSAGES{8} = 'EXIT -- Maximum matrix-vector operations reached';


%----------------------------------------------------------------------
% Initialize iterate history
%----------------------------------------------------------------------

% Pre-allocate iteration info vectors
if (history)
   historyXNorm1 = zeros(min(maxIts,10000),1);
   historyRNorm2 = zeros(min(maxIts,10000),1);
   historyLambda = zeros(min(maxIts,10000),1);
end


%----------------------------------------------------------------------
% Initialization for the hybrid mode
%----------------------------------------------------------------------

% Make sure we can use the hybrid mode
if (hybridMode)
   if ((~realx) || ...
       (~strcmp(func2str(options.project), 'NormL1_project')) || ...
       (~strcmp(func2str(options.primal_norm), 'NormL1_primal')) || ...
       (~strcmp(func2str(options.dual_norm), 'NormL1_dual')))
       Warning('Hybrid mode only applies to non-complex basis pursuit');
       hybridMode = false;
   end
end

% Initialize parameter for the hybrid mode
if (hybridMode)
   % Maintain support
   xAbs   = abs(x);
   xNorm1 = sum(xAbs);
   if ((abs(xNorm1 - tau) / max(1,tau)) < 1e-8)
      support = false(n,1); % Interior
   else
      support  = (xAbs >= 1e-9); % Boundary
   end
   
   H = [];       % Current estimate of the Hessian
   sqrt1 = [];   % Intermediate square-root values #1
   sqrt2 = [];   % Intermediate square-root values #2
   flagUseHessian = false;
   lbfgsHist      = options.lbfgsHist;
else
   flagUseHessian = false;
end


%----------------------------------------------------------------------
% Log header
%----------------------------------------------------------------------
printf('\n');
printf(' %s\n',repmat('=',1,80));
printf(' SPGL1  v.%s (%s)\n', REVISION, DATE);
printf(' %s\n',repmat('=',1,80));
printf(' %-22s: %8i %4s',   'No. rows',       m,      '');
printf(' %-22s: %8i\n',     'No. columns',    n         );
printf(' %-22s: %8.2e %4s', 'Initial tau',    tau,    '');
printf(' %-22s: %8.2e\n',   'Two-norm of b',  bNorm     );
printf(' %-22s: %8.2e %4s', 'Optimality tol', optTol, '');
if singleTau
   printf(' %-22s: %8.2e\n', 'Target one-norm of x', tau);
else
   printf(' %-22s: %8.2e\n', 'Target objective', sigma);
end
printf(' %-22s: %8.2e %4s', 'Basis pursuit tol', bpTol, '');
printf(' %-22s: %8i\n',     'Maximum iterations',maxIts   );
printf('\n');
if singleTau
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %6.1f';
   logH = ' %5s  %13s  %13s  %9s  %6s\n';
   printf(logH,'Iter','Objective','Relative Gap','gNorm','stepG');
else
   logB = ' %5i  %13.7e  %13.7e  %9.2e  %9.3e  %6.1f';
   logH = ' %5s  %13s  %13s  %9s  %9s  %6s  %13s\n';
   printf(logH,'Iter','Objective','Relative Gap','Rel Error','gNorm','stepG','tau');
end


%----------------------------------------------------------------------
% Setup for the first iteration
%----------------------------------------------------------------------

% When sigma >= ||b|| we have x = 0 as the optimal solution. Set tau = 0
% to short-circuit the main loop.
if (~isempty(sigma) && (bNorm <= sigma))
   printf('W: sigma >= ||b||. Exact solution is x = 0.\n');
   tau = 0; singleTau = true;
end

% Project x and compute the residual
x = project(x,tau);
r = b - Aprod(x,1);  % r = b - Ax

% Compute the objective and gradient
f = (r'*r) / 2; 
g = - Aprod(r,2);  % g = -A'r
if (mu > 0)
   f = f + (mu/2) * (x'*x);
   g = g + mu * x;
end

% Required for nonmonotone strategy
lastFv(1) = f;
fBest     = f;
xBest     = x;
fOld      = f;

% Set current maximum dual objective
fDualMax  = -Inf;
gNormBest = NaN;

% Root-finding parameters
tauNew        = tau;
flagFixTau    = false; % Make sure tau is no longer updated
flagUpdateTau = false;

% Compute projected gradient direction and initial steplength.
dx     = project(x - g, tau) - x;
dxNorm = norm(dx,inf); % TODO: Should this be options.norm_dual?
if dxNorm < (1 / stepMax)
   gStep = stepMax;
else
   gStep = min( stepMax, max(stepMin, 1/dxNorm) );
end


%----------------------------------------------------------------------
% MAIN LOOP
%----------------------------------------------------------------------
while 1
    
    %------------------------------------------------------------------
    % Determine the dual objective
    %------------------------------------------------------------------

    gNorm = options.dual_norm(-g,weights); % ||A^Tr - mu*x||
    rtr = (r'*r);
    if (mu > 0), xtx = (x'*x); end

    if (mu == 0)
       % Classic method to determine dual
       fDual = r'*b - tau*gNorm - (rtr)/2;
    else
       % Determine z = |A^Tr| = |-g + mu*x|
       z = abs(mu * x - g);
       
       % Solve the subproblem in lambda:
       %    minimize  tau*lambda + (1/2mu)||[z - lambda*w]_+||_2^2
       objValue = findLambdaStar(z,weightsFull,tau,mu);
       
       % Compute the dual objective
       fDual = r'*b - xtx / 2 - objValue;
    end

    % Use the largest dual objective to determine the relative duality
    % gap for the current primal objective
    if (fDual > fDualMax)
       fDualMax  = fDual;
       gNormBest = gNorm;
    else
       fDual = fDualMax;
    end

    gap = f - fDual;

    % Compute the augmented-residual norm
    if (mu > 0)
       rNorm = sqrt(rtr + mu*xtx);
    else
       rNorm = norm(r, 2);
    end
    
    
    %------------------------------------------------------------------
    % Test exit conditions
    %------------------------------------------------------------------

    % Compute the relative duality gap
    rGap  = abs(gap) / max(relgapMinF,f);
    
    % Too many iterations and not converged.
    if (iter >= maxIts)
       stat = EXIT_ITERATIONS;
    end

    % Check optimality conditions
    if (singleTau)
       if ((rGap <= optTol) || (rNorm < optTol*bNorm))
          stat  = EXIT_OPTIMAL;
       end
    else
       % Test if a least-squares solution has been found -- note that this
       % will be a sub-optimal solution whenever rNorm is less than sigma
       if (gNorm <= lsTol * rNorm)
          stat = EXIT_LEAST_SQUARES;
       end

       % Classic mode
       if (rfModePrimal)
          aError1 = rNorm - sigma;
          aError2 = f - sigma^2 / 2;
          rError1 = abs(aError1) / max(1,rNorm);
          rError2 = abs(aError2) / max(1,f);
          
          if ((rGap <= max(optTol, rError2)) || (rError1 <= optTol))
             if (rNorm <= bpTol * bNorm)
                % Residual minimized: basis-pursuit solution
                stat = EXIT_BPSOL_FOUND;
             elseif (rError1 <= optTol)
                % Found an approximate root
                stat = EXIT_ROOT_FOUND;
             elseif (rNorm <= sigma)
                % Found a sub-optimal basis-pursuit solution
                stat = EXIT_SUBOPTIMAL_BP;
             end
          end
     
          % Check if tau needs to be updated
          testRelChange1 = (abs(f - fOld) <= decTol * f);
          testRelChange2 = (abs(f - fOld) <= 1e-1 * f * (abs(rNorm - sigma)));
          flagUpdateTau  = (((testRelChange1) && (rNorm >  2 * sigma)) || ...
                            ((testRelChange2) && (rNorm <= 2 * sigma))) && ...
                            ~stat && ~flagUpdateTau;
                         
          flagUpdateTau = flagUpdateTau && ~flagFixTau;
                         
          if (flagUpdateTau)
             tauNew = max(0,tau + (rNorm * aError1) / gNormBest); 
          end
       else
          % Dual-based root finding
          
          % When the primal objective is sufficiently close to sigma
          % it is guaranteed that we can solve the problem for the
          % current tau value and we therefore fix tau.
          aError1 = rNorm - sigma;
          rError1 = abs(aError1) / max(relgapMinR,rNorm);
          if (rError1 <= optTol)
             flagFixTau = true;
          end
          
          % The dual objective is used to determine candiate tau values,
          % and we maintain the largest one as the solution for the next
          % root-finding step.
          aErrorDual = (b'*r - tau * gNorm) - rNorm * sigma;
          tauNew = max(tauNew, tau + aErrorDual / gNorm);
          
          if ((flagFixTau) && (rGap <= optTol))
             stat = EXIT_ROOT_FOUND;
          end

          % Root-finding based on the dual
          if ((flagUpdateTau) || (stat) || (flagFixTau))
             flagUpdateTau = false;
          elseif (rGap <= optTol)
             flagUpdateTau = true;
          elseif (((fDual - sigma) / (f - sigma)) > rootfindTol)
             flagUpdateTau = true;
          end
       end
    end

    %------------------------------------------------------------------
    % Update tau if needed
    %------------------------------------------------------------------
    if (flagUpdateTau)
       tauOld   = tau;
       tau      = tauNew;
       nNewton  = nNewton + 1;
       printTau = abs(tauOld - tau) >= 1e-6 * tau; % For log only.
       if (tau < tauOld)
          % The one-norm ball decreased we need to project it to
          % ensure that the next iterate is feasible.
          x = project(x,tau);
             
          % Update the residual, gradient, and function value.
          r = b - Aprod(x,1);  % r = b - Ax
          g =   - Aprod(r,2);  % g = -A'r
          f = r'*r / 2;
          if (mu > 0)
             g = g - mu * x;
             f = f + (mu/2) * (x'*x);
          end

          % Reset the function value history
          lastFv    = -inf(nPrevVals,1);
          lastFv(1) = f;
          fBest     = f;
          xBest     = x;
       end
          
       % Reset Hessian
       if (hybridMode)
          H = []; flagUseHessian = false;
       end
          
       % Reset status
       stat = 0;
          
       % Reset the current maximum dual objective
       fDualMax = -Inf;
       gNormBest = NaN;
    end


    %------------------------------------------------------------------
    % Print log, update history and act on exit conditions.
    %------------------------------------------------------------------
    if ((logLevel >= 2) || singleTau || printTau || (iter == 0) || stat)
       tauFlag = '              ';
       if printTau, tauFlag = sprintf(' %13.7e',tau); end
       if singleTau
          printf(logB,iter,rNorm,rGap,gNorm,log10(stepG));
       else
          printf(logB,iter,rNorm,rGap,rError1,gNorm,log10(stepG));
          if printTau
             printf(' %s', tauFlag);
          end
       end
       printf('\n');
    end
    printTau = false;
    
    % Update history info
    if (history)
       historyXNorm1(iter+1) = options.primal_norm(x,weights);
       historyRNorm2(iter+1) = rNorm;
       historyLambda(iter+1) = gNorm;
    end
    
    if (stat ~= 0), break; end % Act on exit conditions.


    %==================================================================
    % Iterations begin here
    %==================================================================
    iter = iter + 1;
    xOld = x;  fOld = f;  gOld = g;  rOld = r;

    try
       % ===================================
       % Try a quasi-Newton direction first
       % ===================================
       if (flagUseHessian)
          statusMsg = 'H'; iterH = iterH + 1;

          % --------------------------------------------------
          % Step 1. Get search direction
          % --------------------------------------------------
          if (~any(support))
             d = lbfgshprod(H,-g);
          else
             d = -g;

             % Project gradient onto the coefficient space
             dTrans = productBMex(signs .* d(support), 1, sqrt1, sqrt2);

             % If ||dTrans|| is tiny it is (near) orthogonal to the face
             if (norm(dTrans,2) <= 1e-10 * max(1,norm(d,2)))
                %TODO(?): stat = EXIT_OPTIMAL;
                lnErr = true;
             end

             if (~lnErr)
                % Get quasi-Newton search direction
                dQuasi = lbfgshprod(H,dTrans);

                % Convert to global domain to get new search direction
                d = zeros(n,1);
                dSupport =  signs .* productBMex(dQuasi, 0, sqrt1, sqrt2);
                d(support) = dSupport;
             end
          end

          if (~lnErr)
             % --------------------------------------------------
             % Step 2. Determine first non-zero entry to hit zero
             % --------------------------------------------------
             s1 = (((d < 0) & (x > 0)) | ((d > 0) & (x < 0)));
             if (any(s1))
                gamma = -max(x(s1) ./ d(s1));
             else
                gamma = +Inf;
             end

             w = Aprod(d,1);

             %TODO (should be able to omit this)
             lnErr = false;

             % --------------------------------------------------
             % Step 3. Compute the optimal step length beta
             % --------------------------------------------------
             enumerator  = (w'*r);
             denominator = (w'*w);
             if (mu > 0)
                enumerator  = enumerator - mu * (x'*d);
                denominator = denominator + mu * (d'*d);
             end
             beta = enumerator / denominator;
             if (beta <= 1e-11)
                lnErr = true;
             else
                beta = min(beta,gamma);
             end
          end
          
          % --------------------------------------------------
          % Check if a suitable step length was found
          % --------------------------------------------------
          if (~lnErr)
             x = xOld + beta * d;
             r = r - beta * w; % Avoid evaluating A*x
             f = (r'*r) / 2;
             
             if (mu > 0)
                f = f + (mu/2) * (x'*x);
             end

             stepG = beta;
             nLineTot = nLineTot + 1;
             
             % Reset function history (TODO)
             %if (qnResetHist)
                lastFv(:) = f;
             %end
             
             % Record successful QN step (TODO)
             %if (history)
             %   historyQN(iter+1) = true;
             %end
          end

       else
          % Do not use Hessian
          lnErr = true;
       end


       %---------------------------------------------------------------
       % Projected gradient step and linesearch.
       %---------------------------------------------------------------

       % -----------------------------------------------------------------
       % Line search #1: Backtracking curvilinear line search
       % -----------------------------------------------------------------
       if (lnErr)
          [f,x,r,nLine,stepG,lnErr] = ...
              spgLineCurvy(x,gStep*g,max(lastFv),@Aprod,b,@project,tau);
          nLineTot = nLineTot + nLine;

          if (lnErr)
             % Projected backtrack failed. Retry with feasible dir'n linesearch.
             x    = xOld;
             f    = fOld;
             r    = rOld;
          end
       end
       
       % -----------------------------------------------------------------
       % Line search failure: update Barzilai-Borwein scaling factor
       % -----------------------------------------------------------------
       if (lnErr)
          % Revert to previous iterates and damp max BB step.
          x = xOld;
          f = fOld;
          if maxLineErrors <= 0
             stat = EXIT_LINE_ERROR;
          else
             stepMax = stepMax / 10;
             printf(['W: Linesearch failed with error %i. '...
                     'Damping max BB scaling to %6.1e.\n'],lnErr,stepMax);
             maxLineErrors = maxLineErrors - 1;
          end
       end

       % TODO
       ensure(options.primal_norm(x,weights) <= tau+optTol);
       

       %---------------------------------------------------------------
       % Update gradient and compute new Barzilai-Borwein scaling.
       %---------------------------------------------------------------
       if (~lnErr)
          % Compute the gradient
          g = - Aprod(r,2);
          if (mu > 0), g = g + mu * x; end
          
          % Compute the differences
          s    = x - xOld;
          y    = g - gOld;
          sts  = s'*s;
          sty  = s'*y;
          if   sty <= 0,  gStep = stepMax;
          else            gStep = min( stepMax, max(stepMin, sts/sty) );
          end
       else
          gStep = min( stepMax, gStep );
       end
       
       
       %---------------------------------------------------------------
       % Update a Hessian approximation (hybrid mode)
       %---------------------------------------------------------------
       if (hybridMode)
          supportOld = support;
          absx = abs(x); % Used for support determination later
          xNorm1 = sum(absx); % norm(x,1);

          % Step 1. Determine support
          if ((abs(xNorm1 - tau) / max(1,tau)) > 1e-8)
             support = false(n,1);    % Interior of the L1 ball
             nSupport = n+1;
          else
             support  = (absx > 1e-9);
             nSupport = sum(support);
          end
          
          % Step 2. Check if the Hessian should be updated
          if ((nSupport > 1) && (nSupport <= m) && ...
              (iter > 1) && (~lnErr) && ...
              (all(support == supportOld)) && ...
              (all(sign(x(support)) == sign(xOld(support)))))
             % Support remained the same
             flagUpdateHessian = true;
          else
             flagUpdateHessian = false;
          end
          
          % Step 3. Check self-projection condition
          if (flagUpdateHessian)
             if (any(support)) % Boundary of crosspolytope
                d = -g;
                s1 = ((d <  0) & (x > 0)) | ((d >  0) & (x < 0));
                s2 = ((d <= 0) & (x < 0)) | ((d >= 0) & (x > 0));
                s3 = ~(s1 | s2); % Zero entries in x
                
                % Compute the absolute sum of the sets
                sum1 = sum(abs(d(s1))); % Decreases the one norm
                sum2 = sum(abs(d(s2))); % Increases the one norm
                sum3 = sum(abs(d(s3))); % Increases the one norm
                
                if (sum1 > sum2+sum3)
                   % Move into the polytope
                   flagSelfProj = false;
                elseif (sum1 < sum2+sum3)
                   % Move away from the boundary
                   flagSelfProj = (max(abs(d(s3))) <= (sum1+sum2) / sum(support));
                else
                   % Move across the face
                   flagSelfProj = (sum3 == 0);
                end
             else % Interior of crosspolytope
                flagSelfProj = true;
             end
             
             if (~flagSelfProj)
                flagUpdateHessian = false;
             end
          end
          
          % Step 4. Reset or update Hessian approximation
          if (flagUpdateHessian)
             if (isempty(H))
                if (~any(support))
                   % Interior of crosspolytope
                   H = lbfgsinit(n,lbfgsHist,1e-3);
                else
                   % Boundary of crosspolytope
                   signs = sign(x(support));
                   H = lbfgsinit(nSupport-1,lbfgsHist,1e-3);
                end
             end
             
             % Update Hessian approximation
             if (~any(support))
                H = lbfgsupdate(H, 1, s, gOld, g);
             else
                if (isempty(sqrt1))
                   % Compute all square root factors
                   sqrt1 = sqrt(1 ./ ((1:n) .* (2:n+1)))'; % sqrt(1 / (i*(i+1)))
                   sqrt2 = sqrt((1:n) ./ (2:n+1))';        % sqrt(i / (i + 1.0))
                end
                
                H = lbfgsupdate(H, 1, ...
                   productBMex(signs.*s(support)   ,1,sqrt1,sqrt2), ...
                   productBMex(signs.*gOld(support),1,sqrt1,sqrt2), ...
                   productBMex(signs.*g(support)   ,1,sqrt1,sqrt2));
             end
             
             % We can use the Hessian
             flagUseHessian = true;
          else
             flagUseHessian = false;
             H = [];
          end
       end % hybridMode
      
       
    catch err % Detect matrix-vector multiply limit error
       if strcmp(err.identifier,'SPGL1:MaximumMatvec')
         stat = EXIT_MATVEC_LIMIT;
         iter = iter - 1;
         x = xOld;  f = fOld;  g = gOld;  r = rOld;
         break;
       else
         rethrow(err);
       end
    end
    
    
    %------------------------------------------------------------------
    % Update function history.
    %------------------------------------------------------------------
    if ((singleTau) || (f > sigma^2 / 2)) % (Don't update if superoptimal)
       lastFv(mod(iter,nPrevVals)+1) = f;
       if fBest > f
          fBest = f;
          xBest = x;
       end
    end

end % while 1

%----------------------------------------------------------------------
% END OF MAIN LOOP
%----------------------------------------------------------------------

% Restore best solution (only if solving single problem).
if ((singleTau) && (f > fBest))
   rNorm = sqrt(2*fBest);
   printf('\n Restoring best iterate to objective %13.7e\n',rNorm);
   x = xBest;
   r = b - Aprod(x,1);
   g =   - Aprod(r,2);
   if (mu > 0)
      g = g + mu * x;
      rNorm = sqrt(r'*r + mu * x'*x);
   else
      rNorm = norm(r,2);
   end
   gNorm = options.dual_norm(g,weights);
end

% Final cleanup before exit.
info.tau         = tau;
info.rNorm       = rNorm;
info.rGap        = rGap;
info.gNorm       = gNorm;
info.rGap        = rGap;
info.stat        = stat;
info.iter        = iter;
info.nProdA      = nProdA;
info.nProdAt     = nProdAt;
info.nNewton     = nNewton;
info.timeProject = timeProject;
info.timeMatProd = timeMatProd;
info.options     = options;
info.timeTotal   = toc;

if (history)
   info.xNorm1      = historyXNorm1(1:iter);
   info.rNorm2      = historyRNorm2(1:iter);
   info.lambda      = historyLambda(1:iter);
end

% Print final output
if ((stat < 1) && (stat > length(EXIT_MESSAGES)))
   error('Unknown termination condition\n');
end
printf('\n %s\n', EXIT_MESSAGES{stat})
printf('\n');
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A',nProdA,'','Total time   (secs)', info.timeTotal);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Products with A''',nProdAt,'','Project time (secs)',timeProject);
printf(' %-20s:  %6i %6s %-20s:  %6.1f\n',...
   'Newton iterations',nNewton,'','Mat-vec time (secs)',timeMatProd);
printf(' %-20s:  %6i\n', ...
   'Line search its',nLineTot);
printf('\n');



% =====================================================================
% NESTED FUNCTIONS.  These share some vars with workspace above.
% =====================================================================
    
function z = Aprod(x,mode)
   if (nProdA + nProdAt >= maxMatvec)
     error('SPGL1:MaximumMatvec','');
   end
     
   tStart = toc;
   if mode == 1
      nProdA = nProdA + 1;
      if   explicit, z = A*x;
      else           z = A(x,1);
      end
   elseif mode == 2
      nProdAt = nProdAt + 1;
      if   explicit, z = A'*x;
      else           z = A(x,2);
      end
   else
      error('Wrong mode!');
   end
   timeMatProd = timeMatProd + (toc - tStart);
end % function Aprod

% ----------------------------------------------------------------------
function printf(varargin)
% ----------------------------------------------------------------------
  if logLevel > 0
     fprintf(fid,varargin{:});
  end
end % function printf


% ----------------------------------------------------------------------
function x = project(x, tau)
% ----------------------------------------------------------------------
   tStart = toc;
   x = options.project(x,weights,tau);
   timeProject = timeProject + (toc - tStart);
end % function project

% =====================================================================
% End of nested functions.
% =====================================================================

end % function spg


% =====================================================================
% PRIVATE FUNCTIONS
% =====================================================================

% ----------------------------------------------------------------------
function [fNew,xNew,rNew,iter,err] = spgLine(f,x,d,gtd,fMax,Aprod,b)
% ----------------------------------------------------------------------
% Nonmonotone linesearch.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
maxIts = 10;
step   = 1;
iter   = 0;
gamma  = 1e-4;
gtd    = -real(gtd);

while 1

    % Evaluate trial point and function value.
    xNew = x + step*d;
    rNew = b - Aprod(xNew,1);
    fNew = rNew'*rNew / 2;

    % Check exit conditions.
    if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif  iter >= maxIts           % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    
    % Safeguarded quadratic interpolation.
    if step <= 0.1
       step  = step / 2;
    else
       tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
       if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
          tmp = step / 2;
       end
       step = tmp;
    end
    
end % while 1

end % function spgLine


% ----------------------------------------------------------------------
function [fNew,xNew,rNew,iter,step,err] = ...
    spgLineCurvy(x,g,fMax,Aprod,b,project,tau)
% ----------------------------------------------------------------------
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
EXIT_NODESCENT  = 2;
gamma  = 1e-4;
maxIts = 10;
step   =  1;
sNorm  =  0;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
iter   =  0;
n      =  length(x);
   
while 1

    % Evaluate trial point and function value.
    xNew     = project(x - step*scale*g, tau);
    rNew     = b - Aprod(xNew,1);
    fNew     = rNew'*rNew / 2;
    s        = xNew - x;
    gts      = scale * real(g' * s);
    if gts >= 0
       err = EXIT_NODESCENT;
       break
    end
    
    if fNew < fMax + gamma*step*gts
       err = EXIT_CONVERGED;
       break
    elseif iter >= maxIts % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search
    % directions can give exactly the same point after projection.  If
    % we observe this in adjacent iterations, we drastically damp the
    % next search direction.
    % 31 May 07: Damp consecutive safeguarding steps.
    sNormOld  = sNorm;
    sNorm     = norm(s) / sqrt(n);
    %   if sNorm >= sNormOld
    if abs(sNorm - sNormOld) <= 1e-6 * sNorm
       gNorm = norm(g) / sqrt(n);
       scale = sNorm / gNorm / (2^nSafe);
       nSafe = nSafe + 1;
    end
    
end % while 1

end % function spgLineCurvy
