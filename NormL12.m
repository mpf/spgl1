classdef NormL12 < NormObj
   properties
      blocksize = 1
   end
   
   methods
      function obj = NormL12(blocksize, weights)
         obj = obj@NormObj();
         obj.blocksize = blocksize;
         
         if (nargin == 2)
            obj.weights = weights;
         end
      end
      
      function p = primal(obj, x)
          m = round(length(x) / obj.blocksize); n = obj.blocksize;
          if isreal(x)
             p = sum(obj.weights.*sqrt(sum(reshape(x,m,n).^2,2)));
          else
             p = sum(obj.weights.*sqrt(sum(abs(reshape(x,m,n)).^2,2)));
          end
      end
      
      function d = dual(obj, y)
         m = round(length(y) / obj.blocksize); n = obj.blocksize;
         if isreal(y)
            d = norm(sqrt(sum(reshape(y,m,n).^2,2))./obj.weights,inf);
         else
            d = norm(sqrt(sum(abs(reshape(y,m,n)).^2,2))./obj.weights,inf);
         end
      end
      
      function p = project(obj, x, tau)
         % Convert to matrix
         m = round(length(x) / obj.blocksize); n = obj.blocksize;
         x = reshape(x,m,n);

         % Compute two-norms of rows
         if isreal(x)
            xa  = sqrt(sum(x.^2,2));
         else
            xa  = sqrt(sum(abs(x).^2,2));
         end

         % Project one one-norm ball
         idx = xa < eps;
         xc  = oneProjector(xa,obj.weights,tau);

         % Scale original
         xc  = xc ./ xa; xc(idx) = 0;
         p   = spdiags(xc,0,m,m)*x;

         % Vectorize result
         p = p(:);
      end
   end
end
