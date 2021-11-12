function bool = inrange(x,range)
% inrange(x,range) returns True if x is in the interval range

  % expand the range a little for rounding errors...
  range = sort(range) + 1e4*[-eps,eps];
  bool = logical(mod(sum(x<=range),2));
  return
  