% ### RELAXdynamic.m ###       04.29.10

% simple code to help figure out how to design a stimulus waveform for
% essentially tone bursts (with some steady-state region) to determine
% effects upon SOAEs (would need to entirely re-code this in C)

% this version uses an exponential to get the desired sigmoidal behavior in
% the envelope (coupled with a conditional statement re the specified
% times); intial testing indicates it seems to work well

% NOTE: this could likely be modified without too much difficulty to
% simultaneously present two tones


clear

% ------------------------

fSTIM= 1200;   % [kHz]
LevSTIM= 1;    % stim. aplitude [arb. at this point]

tauON= 0.5;   % [ms] time to come to gamma% (see below) of the full amplitude
tauOFF= 0.5;   % [ms] time to go from full amplitude to 1-gamma% 

delayON= 10;    % [ms] delay for ramp ON
duration= 50;   % [ms] duration signal is ON 'full' before OFF envelope starts
delayOFF= 50;   % [ms] after signal is OFF, duration of zero padding

repeatN= 3;    % # of times to repeat

gamma= 0.9;   % fraction of amplitude signal should be at by tauON

phaseON= 0;   % phase of sinusoid at onset [-pi,pi] NOTE: code uses sin, not cos

SR= 44100;  % [Hz]

Npoints= 8192;     % length of fft window (# of points) [should ideally be 2^N]
                   % [not currently used, but potentially helpful for
                   % post-analysis]

% -----------------------


% ^^^^^^^^^^^^^^^^^^^
% create array of time values
lengthT= delayON+tauON+duration+tauOFF+delayOFF;  % total window length [ms]
t= [0:1/SR:lengthT/1000];


% ^^^^^^^^^^^^^^^^^^^
dt= 1/SR;  % spacing of time steps
% create a freq. array (for FFT bin labeling)
freq= [0:Npoints/2];
freq= SR*freq./Npoints;



% ^^^^^^^^^^^^^^^^^^^
% quantize the freq. (so to have an integral # of cycles in FFT window)
df = SR/Npoints;
fQ= ceil(fSTIM/df)*df;   % quantized natural freq.

disp(sprintf('specified freq. = %g Hz', fSTIM));
disp(sprintf('quantized freq. = %g Hz', fQ));
fSTIM= fQ;  % reset feq. as quantized version


% ^^^^^^^^^^^^^^^^^^^
% determine carrier signal (such that phase is zero at envelope onset,
% uless specified otherwise)
carrier= sin(2*pi*fSTIM*(t-delayON/1000) - phaseON);


% ^^^^^^^^^^^^^^^^^^^
% convert specified times to s
delayON= delayON/1000;
tauOFF= tauOFF/1000;
tauON= tauON/1000;
duration= duration/1000;


% ^^^^^^^^^^^^^^^^^^^
% determine constant for onset & offset envelopes
a= -log(1-gamma)/tauON;
b= -log(1-gamma)/tauOFF;


% ^^^^^^^^^^^^^^^^^^^
% determine envelope
for nn=1:size(t,2)
    if t(nn)<delayON
        envelope(nn)= 0;
    % determine onset envelope + 'on' envelope  
    elseif t(nn)>= delayON && t(nn)<delayON+tauON+duration/2
        envelope(nn)= LevSTIM*(1-exp(-a*(t(nn)-delayON)));
    % determine offset envelope (NOTE: think about SYMMETRY to get this correct)
    elseif t(nn)>= delayON+tauON+duration/2 && t(nn)< delayON+tauON+duration+tauOFF
        %%envelope(nn)= LevSTIM*(1-exp(-b*(t(nn)-delayON-tauON-duration/2)));
        envelope(nn)= LevSTIM*(1-exp(b*(t(nn)-delayON-tauON-duration-tauOFF)));
    elseif t(nn)>= delayON+tauON+duration+tauOFF
        envelope(nn)= 0;
    end
end



% ^^^^^^^^^^^^^^^^^^^
% apply onset envelope to carrier
% [NOTE: this is only for the initial window only]
signal= envelope.*carrier;


% ^^^^^^^^^^^^^^^^^^^
% copy/paste as needed to produce the repeated segements
signalOG= signal;           % set aside the original signal
indx= size(signal,2)+1;     % required indexing offset
offset= size(t,2);
for nn=1:(repeatN-1)*offset
    t(indx)= t(indx-1)+dt; 
    signal(indx)= signal(indx-offset);
    indx=indx+1;    %increment counter
end


% ^^^^^^^^^^^^^^^^^^^
% determine spectra of original signal
spec= rfft(signalOG);

Npoints2= size(signalOG,2);  % # of time points in waveform
% determine corresponding freqs. for FFT bins
freq= [0:Npoints2/2];
freq= SR*freq./Npoints2;

figure(2); clf;
plot(freq/1000,db(spec),'kx-')
grid on; hold on;
axis([0 6 -120 0])
xlabel('freq. [kHz]')
ylabel('amplitude [dB]')
title('Spectrum of tone w/ envelope (single presentation)')


% ^^^^^^^^^^^^^^^^^^^
% visualize the time waveform
figure(1); clf;
%plot(t,carrier);
hold on; grid on;
plot(t(1:offset),envelope,'r.-')
plot(t(1:offset),-envelope,'r-')
plot(t,signal,'b.-');
xlabel('time [ms]')
ylabel('output [e.g., V]')


% *****************************************
% play the stimuli
if 1==0
    sound(signal,SR)
end
