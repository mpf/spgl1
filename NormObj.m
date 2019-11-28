classdef (Abstract) NormObj
   properties
      weights = 1;
   end
   
   methods
      function p = primal(obj,x)
      end
      
      function d = dual(obj,y)
      end

      function p = project(obj,x, tau)
      end
   end
end
