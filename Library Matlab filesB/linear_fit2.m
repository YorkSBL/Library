function [a,b,Sigma_a,Sigma_b] = linear_fit2 (x,y,Sigma_y)
%
% [a,b,Sigma_a,Sigma_b,Q,Cov_ab,R_ab,Chi2,R] = linear_fit (x,y,?Sigma_y?)
%
% Given data vectors (x,y) with optional individual standard 
% deviations, Sigma_y, fits them to a straight line, 
%         y = a + b*x.
% Returns intercept a and slope b and their respective probable
% uncertainties, Sigma_a and Sigma_b. 
%
% [modified/simplified version of linear_fit.m as written by C. Shera;
% main difference here is that I am using formulae from Bevington that I 
% am a bit more familiar with (should lead to the same #s!!)
%
% C. Bergevin
%
estimate_errors = 0;
if (nargin < 3)
  Sigma_y = ones (size (y));
  estimate_errors = 1;
end
if (length (x) ~= length (y))
  error ('Vectors x and y must have the same length');
end
if (length (x) < 2)
  error ('Need at least two points to fit a line');
end
if (length (Sigma_y) == 1)
  Sigma_y = Sigma_y * ones (size (y));
end
if (length (Sigma_y) ~= length (y))
  error ('Vectors y and Sigma_y must have the same length');
end

x       = x(:);
y       = y(:);
Sigma_y = Sigma_y(:);

wt = 1./Sigma_y.^2;

S  = sum (wt);
Sx = sum (x.*wt);
Sy = sum (y.*wt);

Sxx = sum (x.*x.*wt);
Sxy = sum (x.*y.*wt);

Delta = S*Sxx - Sx^2;
coeff = 1/Delta;

N        = length (x);

b = coeff*(S*Sxy - Sx*Sy);
a = coeff*(Sxx*Sy - Sx*Sxy);


Sigma_a  = sqrt (coeff*Sxx);
Sigma_b  = sqrt (coeff*S);


% if (estimate_errors)
%   Q = 1.0;		  % assume good fit
%   factor = sqrt(Chi2/(N-2));
%   Sigma_a  = Sigma_a * factor;
%   Sigma_b  = Sigma_b * factor;
% end

if (nargout==1)
  % return everything as a row vector 
  a = [a,b,Sigma_a,Sigma_b];
end


