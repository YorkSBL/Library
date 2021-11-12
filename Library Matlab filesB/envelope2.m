

% ENVELOPE  ( MatLinks) Finds the envelope of an arbitrary function.
%
%    PP=ENVELOPE(Y) returns the pp-form of the cubic spline interpolants of the
%    envelope, for later use with PPVAL, etc.  PP will have two columns, the
%    first of which is the upper envelope of Y, and the second of which is the
%    lower envelope.
%
%    YY=ENVELOPE(Y,XX) interpolates (both the upper and lower) envelope of Y at
%    the points specified by abscissa XX, and returns the envelope as the vector
%    YY.  YY has two columns, the first of which is the upper envelope of Y, and
%    the second of which is the lower envelope.
%
%    ENVELOPE(...) by itself plots Y and its envelope of oscillation.
%
%    See also PPVAL, BEATFREQ.
%
%    Type HELP MATLINKS for a full listing of all MatLinks ToolChest functions.
%
function [yy] = envelope(y, xx)
%===============================================================================
%  Copyright  1998,2000 Julian Andrew de Marchi, Ph.D. (julian@matlinks.net)
%  Use & distribution covered by GNU General Public License (www.gnu.org)
%===============================================================================

%------------------
% parse the inputs
%------------------
if (nargin == 0), error('No data vector X supplied.'); end;

%-----------------------
% measure the envelopes
%-----------------------
[i,j] = findpeak(y);

% functional form
if (nargin == 1),

  yy = [spline(i, y(i))' spline(j, y(j))'];
  xx = merge(i,j);
  y1 = ppval(yy(:,1), xx);
  y2 = ppval(yy(:,2), xx);

% vector form
else

  y1 = spline(i, y(i), 1:length(xx))';
  y2 = spline(j, y(j), 1:length(xx))';

  if (length(y1) > length(y2)),
    y2 = [y2; zeros(length(y2) - length(y1), 1)];
  elseif (length(y2) > length(y1)),
    y1 = [y1; zeros(length(y1) - length(y2), 1)];
  end;

  yy = [y1 y2];

end;


%--------------------------------------------------
% plot the envelopes if there's no output variable
%--------------------------------------------------
if (nargout==0),

  hold off, plot(xx, interp1(y, 1:length(xx), '*linear')),
  hold on, plot(xx,y1,'g:', xx,y2,'c:'),
  xlabel('i'), ylabel('y(i)'), title('Envelope'), zoom on;

end;

%===============================================================================
% End-of-File
%===============================================================================

