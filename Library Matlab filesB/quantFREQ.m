function y = quantFREQ (f,SR,N)
% provides the correct quantized freq. such that an integral # of periods
% will fit into an FFT window
% C. Bergevin

df= SR/N;
y= ceil(f/df)*df;


return




