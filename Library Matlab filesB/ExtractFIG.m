
% 09.29.11

% Simple code that demonstrates the syntax to extract data 
% fields from a messy .fig file; note you can also use CS' 
% extract.m to semi-manually do this

clear

file= 'coch_v1_Q20_onset2ms_fFixed1000Hz.fig';  % specify .fig file

hndl= openfig(file);    % opens .fig file (and loads in data from memory)


pH= subplot(414)    % grabs the bottom phase plot, set pH as the handle
grid on;

% -------
% ** to do this manually **
% 1. click on the curve you want
% 2. type > [y,x]=extract;
% 3. now you can plot, e.g., > plot(x,y,'b')

% -------
% ** to do this (semi-kludge) automatically **
t1= findall(pH,'Type','line');      % creates new handle set specific to x/y data

% now the commands below can extract out the desired data 
% 1. grab CF location 
x1= get(t1(1),'xdata');
y1= get(t1(1),'ydata');
% 2. grab peak location 
x2= get(t1(2),'xdata');
y2= get(t1(2),'ydata');
% 1. grab phase data 
x3= get(t1(3),'xdata')
y3= get(t1(3),'ydata')


figure(4); clf;
plot(x3,y3);
hold on; grid on;


