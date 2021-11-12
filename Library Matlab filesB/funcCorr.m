function [Data] = funcCorr(P);
% *** funcCorr.m ***       2018.03.19       C. Bergevin

% [IN PROGRESS; incomplete]
% Function to correlate two (1-D) waveforms in a few different ways

% In.type=1 - simple straight-up correlation
% In.type=2 - create a periodic boundary condition (e.g., via repmat) for
% one of them so to avoid roll-off when multiplying zero elements
% In.type=3 - similar to type=2 but instead do correlation over a
% (user-specified) smaller region of one of the waveforms

% Notes
% o there is a CB-crafted correlation routine (correlate1.m) which is used
% in several example codes (e.g., EXautocorr5.m, EXhoResonanceC.m), but is
% relatively slow compared to Matlab's built-in xcorr.m, so will instead
% aim to use that latter (blackbox) code