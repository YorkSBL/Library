function y = makeSpectrogram(file)

% ### makeSpectrogram.m ###         11.13.12

% reads in data created via RECspeech.cpp and makes a spectrogram

% NOTE: make sure sample rate specified here matches that used when
% recording the data!

SR= 44100;  % SR data collected at [Hz]
disp(sprintf('assumed sample rate = %g kHz', SR/1000));

A= load(file);

%plot(A)
spectrogram(A,128,0,1024,SR,'yaxis');
axis([0 size(A,1)/SR 0 10000])
colorbar

% save data to .wav file?
if 1==0
    fileN= [file(1:end-3),'wav'];
    wavwrite(A,SR,16,fileN);
end

% NOTE: to play back the .wav file (via Matlab), type:
% > B=wavread('AKname1.wav');
% > sound(B,SR)
% where SR is the appropriate sample rate (e.g., fiddle with if you want to
% change the pitch)