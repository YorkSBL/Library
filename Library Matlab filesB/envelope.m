function y = envelope(x)
%ENVELOPE   Calculates signal envelope
%   y = envelope(x) calculates envelope of the input signal x using Hilbert
%   transform (in Signal Processing Toolbox).
%
%   If you want to show an envelope shape of a real-world audio signals,
%   you might want to consider applying some low-pass filter on the
%   envelope to show smooth shape.
%
%   Example:
%      fs = 44100;
%      t = 0:1/fs:0.02;
%      x = sin(2*pi*100*t) + 1/4 * sin(2*pi*1000*t); % 100 & 1000Hz signals
%      y = envelope(x);
%      plot(t,x,'c', t,y,'k') % superippose two signals
%
%   2002-10-?? by MARUI Atsushi
%   2004-10-30 added some documentation

y = abs(hilbert(x));

