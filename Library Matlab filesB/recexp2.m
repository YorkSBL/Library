function Sn = recexp2(f,CF,BW,n)
% modified from recexp.m (C. Shera) so to allow filter
% to act in freq. domain, given a CF and BW


% Recursive exponential filtering.
% Arguments:
%   f  -- "frequency" axis (array)
%   CF -- center "frequency" (scalar)
%   BW -- filter bandwidth (scalar)
%   n  -- filter order
% See Shera and Zweig (1993) J. Acoust. Soc. Am, 93:3333-3352.
% Appendix Eqs. A4-A6 

% Example:
%   freq=linspace(0,11050,4092);
%   plot(freq,recexp2(freq,5000,500,10));  
% 
% Christopher A. Shera (shera at mit.edu)
% Eaton-Peabody Laboratory  
  
  if (nargin<4), n=10; end;

  gamma_n = 1;
  % compute lambda_n...
  for k=2:n
    gamma_n = log(gamma_n+1);
  end
  lambda_n = sqrt(gamma_n);
    
  f= f-CF;  % shift filter to specified CF
  
  f = lambda_n*f/BW;

  Gamma_n = exp(f.^2);
  % compute Gamma_n
  for k=2:n
    Gamma_n = exp(Gamma_n-1);
  end
  
  Sn = 1./Gamma_n;
  return
  
  