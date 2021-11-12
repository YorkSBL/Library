function [delta,rho,mu] = dpole (y);
%
% function [delta,rho,mu] = dpole (y);
%     or
% function params = dpole (y);
%     returns matrix whose columns are values of delta,rho,mu
% function dpole (y);
%     writes to global variables delta,rho,mu
%
% Compute parameter values for Zweig active model
% with a double pole near beta=1 at a distance y 
% from the real axis. See S&Z 1995 footnote 8

if (nargout == 0)
  global delta rho mu;
end

% we hard-wire these for speed
% epsilon for n=2 gives mu near 1.75
% (see epsilon.m for its computation)
% epsilon = 0.09170844119620;
% c = 2 + 1/(epsilon*epsilon)
c = 120.8998691636393;

if (isrow(y))
  y = y(:);
end

yy = y.*y;

a = (y + sqrt (yy + c*(1-yy))) / c;
delta = 2*(y-a);
mu = 1./(2*pi*a);
pimu = 1./(2*a);                     % = pi*mu;
rho = sqrt (1 - delta.*delta/4) ./ (pimu.*exp (1+pimu.*delta));

if (nargout == 1)
   delta = [delta,rho,mu];
end
