function  [err, J] = Dapprox(fun,x,d,f,varargin)
% Approximate Jacobian by forward differences

% This file is part of the "immoptibox" package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package.

J = zeros(length(f),length(x));
xx = x;
for  j = 1 : length(x)
  if  x(j) == 0,  xp = d^2;
  else,           xp = x(j) + d*abs(x(j)); end
  xx(j) = xp;  fp = feval(fun,xx,varargin{:});
  J(:,j) = (fp - f)/(xp - x(j));
  xx(j) = x(j);
end
% Check J
if  ~isreal(J) | any(isnan(J(:))) | any(isinf(J(:)))
  err = -6;  
else,  err = 0;  end