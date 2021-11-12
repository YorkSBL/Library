function x = logscale (x0, x1, npts, rowvec)
% function x = logscale (x0, x1, ?npts=128?,?rowvec=true?);)
% Create logarithmically spaced vector with npts elements ranging from x0 to x1. 
% What logspace should have done.
%
% Christopher A. Shera (shera at mit.edu)
% Eaton-Peabody Laboratory  
%

if (nargin < 3), npts = 128; end;
if (nargin < 4), rowvec = true; end;

x = logspace (log10(x0),log10(x1),npts);

% make sure the end points match exactly
% [They are often off by eps() or so, 
% which can cause havoc for interp1 and other programs.]
x(1) = x0;
x(end) = x1;

if (~rowvec)
  x = x(:);
end
