function y= findpeaksCB(A,In);        % 08.09.15  (C. Bergevin)

% custom version of code to "find peaks" (initially in context of SOAE
% specta) to improve upon the deficiencies/limitations of Matlab's built-in
% function findpeaks.m

% NOTES:
% o If In.thresh<0, this can be used to find "dips/valleys"
% o Not terribly good when two values straddle either side of a peak... :-(

% Input
% A: 1-D waveform (Nx1 array)
% In: structure containing various options (see below)


if (~isfield(In,'thresh')), In.thresh= 1; end  % # of units that a point must have about neighbors to be considered a peak/dip
if (~isfield(In,'range')), In.range= 3; end  % # of (plus or minus index) neighbors to include for comparison
if (~isfield(In,'smooth')), In.smooth= 0; end  % (boolean) use loess to locally smooth (in case points are "jumpy")?


store= [];
for nn= In.range+1:numel(A)-In.range-1
    
    % build up logical arrays to check that all points "passed"
    cc= (A(nn)-A(nn-In.range:nn-1)>=In.thresh); % lower indicies
    dd= (A(nn)-A(nn+1:nn+In.range)>=In.thresh); % higher indicies
    temp= [cc; dd];
    % store away indicies that passed (uses Matlab's function all.m)
    if (all(temp)), store= [store nn];  end
    
end

y= store;

return