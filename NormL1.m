classdef NormL1 < NormObj
   methods
      function obj = NormL1(weights)
         obj = obj@NormObj();
         
         if (nargin == 1)
            obj.weights = weights;
         end
      end
      
      function p = primal(obj, x)
         p = norm(x.*obj.weights,1);
      end
      
      function d = dual(obj, y)
         d = norm(y./obj.weights,inf);
      end
      
      function p = project(obj, x, tau) 
         if isreal(x)
            p = oneProjector(x,obj.weights,tau);
         else
            xa  = abs(x);
            idx = xa < eps;
            xc  = oneProjector(xa,obj.weights,tau);
            xc  = xc ./ xa; xc(idx) = 0;
            p   = x .* xc;
         end
      end
   end
end
