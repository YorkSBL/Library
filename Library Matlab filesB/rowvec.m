function rv = rowvec(v)
% function rv = rowvec(v)
% Convert v to a rowvector
  
  rv = reshape(v,1,prod(size(v)));
  return
  