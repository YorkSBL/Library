function [y,n] = ininterval (x,interval,flag)
%ininterval: returns members of x falling within a specified interval 
% function [y,n] = ininterval (x,interval,?flag=0?)
% By default (flag=0), returns those elements of the vector x that
% fall within the interval (e.g, y=ininterval(11:20,[12,14]) returns 12,13,14.
% If flag==1, returns the indices rather than the elements themselves.
% If flag==2, returns a boolean vector with the same length as x
% If nargout>1, returns in n the indices
  
  if (nargin<3), flag = 0; end;
  
  % expand the interval a little for rounding errors...
  interval = interval + 1e4*[-eps,eps];
  bool = logical(x>=interval(1) & x<=interval(2));

  n = find(bool);
  if (flag == 0)
    y = x(x>=interval(1) & x<=interval(2));
  elseif (flag == 1)
    y = n;
  elseif (flag ==2)
    y = bool;
  end
  
  



