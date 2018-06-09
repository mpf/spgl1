function varargout = findLambdaStar(z,w,tau,mu)
% Solves:
%    minimize  tau*lambda + (1/(2*mu)) * ||[z - lambda*w]_+ ||_2^2
%     lambda
% Returns:
%    Optimal lambda, objective value

lambdak = (z ./ w)';
[lambdak,idx] = sort(lambdak,'descend');
w = w(idx);
z = z(idx);

g  = tau;
i  = 0;
ww = 0;
wz = 0;
n  = length(lambdak);

while (g > 0) && (i < n)
   i = i + 1;
   ww = ww + w(i)*w(i);
   wz = wz + w(i)*z(i);
   
   if (i < n)
      lambda = lambdak(i+1);
   else
      lambda = 0;
   end
   
   g = tau - (wz - lambda * ww) / mu;
end

% Locally we have gradient tau - (1/mu)*(wz - lambda * ww))
% solve for lambda
lambdaStar = max(0, wz - tau*mu) / ww;

c   = max(z - w*lambdaStar,0);
obj = tau * lambdaStar + (1/(2*mu)) * (c'*c);

if (nargout == 1)
   varargout{1} = obj;
else
   varargout{1} = lambdaStar;
   varargout{2} = obj;
end
