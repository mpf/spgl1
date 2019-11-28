classdef NormGroupL2 < NormObj
   properties
      groups = [];
   end
      
   methods
      function obj = NormGroupL2(groups, weights)
         obj = obj@NormObj();
         
         g = groups(:);
         n = length(g);
         [gidx,idx1,idx2] = unique(g);
         obj.groups = sparse(idx2,1:n,ones(1,n),length(gidx),n);
         
         if (nargin == 2)
            obj.weights = weights;
         end
      end
      
      function p = primal(obj, x)
         if isreal(x)
            p = sum(obj.weights.*sqrt(sum(obj.groups * x.^2,2)));
         else
            p = sum(obj.weights.*sqrt(sum(obj.groups * abs(x).^2,2)));
         end
      end
      
      function d = dual(obj, y)
         if isreal(y)
            d = norm(sqrt(sum(obj.groups * y.^2,2))./obj.weights,inf);
         else
            d = norm(sqrt(sum(obj.groups * abs(y).^2,2))./obj.weights,inf);
         end
      end
      
      function p = project(obj, x, tau)
         % Compute two-norms of rows
         if isreal(x)
            xa = sqrt(sum(obj.groups * x.^2,2));
         else
            xa = sqrt(sum(obj.groups * abs(x).^2,2));
         end

         % Project one one-norm ball
         idx = xa < eps;
         xc  = oneProjector(xa,obj.weights,tau);

         % Scale original
         xc  = xc ./ xa; xc(idx) = 0;
         p   = full(obj.groups' * xc).*x;
      end
   end
end
