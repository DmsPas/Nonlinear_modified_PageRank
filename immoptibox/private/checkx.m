function  [err, x,n] = checkx(x0)
%CHECKX  Check vector

% Version 04.01.25.

% This file is part of the "immoptibox" package and is distributed under
% the 3-Clause BSD License. A separate license file should be provided as
% part of the package.

err = 0;  sx = size(x0);   n = max(sx);
if  (min(sx) ~= 1) | ~isreal(x0) | any(isnan(x0(:))) | isinf(norm(x0(:))) 
  err = -1;   x = []; 
else
  x = x0(:); 
end