function Sn = recexp (t,tc,n)
% function Sn = recexp (t,tc,?n=10?)
% Recursive exponential filtering.
% Arguments:
%   t  -- "time" axis (array)
%   tc -- cutoff "time" (scalar)
%   n  -- filter order
% See Shera and Zweig (1993) J. Acoust. Soc. Am, 93:3333-3352.
% Appendix Eqs. A4-A6
% Example:
%   t = linspace(0,1,256);
%   plot(t,recexp(t,0.3,10));  
% 
% Christopher A. Shera (shera at mit.edu)
% Eaton-Peabody Laboratory  
  
  if (nargin<3), n=10; end;

  gamma_n = 1;
  % compute lambda_n...
  for k=2:n
    gamma_n = log(gamma_n+1);
  end
  lambda_n = sqrt(gamma_n);

  t = lambda_n*t/tc;

  Gamma_n = exp(t.^2);
  % compute Gamma_n
  for k=2:n
    Gamma_n = exp(Gamma_n-1);
  end
  
  Sn = 1./Gamma_n;
  return
  
  