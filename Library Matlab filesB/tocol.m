function y = tocol(x)
% function y = tocol(x)
% return x as a column vector
  
  [m,n] = size(x);
  if (m ~= 1 & n ~= 1)
    error ('x not a vector');
  end

  y = x(:);
  
  return
  