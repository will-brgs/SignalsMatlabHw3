%% Signals and Systems Matlab Homework #3
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 351
% * Date:                     Created 2/11/2024, Last Edited 2/20/2024
% * With contributions from:  Mack Larosa, Tasha Igic, Mischa Tranor
%% Initialize Variables
R = 1e3; %Resistance in ohms
C = 5e-6; %Capacitence in Farads
tau = R*C;
sampleFreq = 44.1e3;
samplePeriod = 1/sampleFreq;

%% Part 1: Complex Exponential Input
freqVals = [10, 1000];

for i = 1:2
freq = freqVals(i);
anglularFreq = 2*pi*freq;

inputFunct = zeros(1,15*tau);
inputFunct = 

end
