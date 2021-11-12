function xs = xings (y,x,flag,tol)
% function xs = xings (y,?x?,?flag='+'?,?tol=20*eps?)
% Extract values of x at zero crossings of the vector y.
% Uses linear interpolation. If isempty(x), uses point number.
% If flag is '+', returns positive-going zero crossings.
% If flag is '-', returns negative-going zero crossings.

if (nargin<2 | isempty(x)), x=(1:length(y)); end;
if (nargin<3 | isempty(flag)), flag='+'; end;
if (nargin<4), tol = 20*eps; end;
if (isstr(flag) & strcmp(flag,'-'))	% negative going
  y = -y;
end
y = y(:);
x = x(:);

ndx = find ((y((1:end-1)')<=0) & (y((2:end)')>0));
xs = zeros(size(ndx));			% preallocate

for i = 1:length(ndx)
  r = (ndx(i):ndx(i)+1)';
  xs(i) = interp1q (y(r),x(r),0);
end

% deal with the ends...
% if the value is close enough to zero (and wouldn't
% already have been detected) we call it a zero...
if (y(1)>0 & y(1)<tol)
  xs = [x(1);xs];
end
if (y(end)<0 & abs(y(end))<tol)
  xs = [xs;x(end)];
end


