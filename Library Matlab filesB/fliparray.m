function x=fliparray(x)
%FLIPARRAY Reverse Array Element Order.
% FLIPARRAY(X) returns an array the same size as X with its elements in
% reverse order. All data types are supported.
%
% FLIPARRAY(X) is a generalization of FLIPLR and FLIPUD.
% When X is a vector it reverses the order of elements no matter if it is a
% row or column.
% When X is a general array, it simply reverses the order of the elements.
%
% See also FLIPLR, FLIPUD, FLIPDIM, ROT90.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2005-08-25

x(:)=x(end:-1:1); % simple one-liner, but useful.