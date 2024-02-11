%% Signals and Systems Matlab Homework #3
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 351
% * Date:                     Created 2/11/2024, Last Edited 2/20/2024
% * With contributions from:  Mack Larosa, Tasha Igic, Mischa Tranor
%% Initialize Variables
close all
R = 1e3; %Resistance in ohms
C = 5e-6; %Capacitence in Farads
tau = R*C;
sampleFreq = 44.1e3;
samplePeriod = 1/sampleFreq;

%% Part 1: Complex Exponential Input
freqVals = [10, 1000];

for i = 1:2
freq = freqVals(i);
angularFreq = 2*pi*freq;

sampleTimes = 0:samplePeriod:15*tau;
%inputFunct = zeros(1,length(timePoints));
inputFunct = exp(1j*angularFreq*sampleTimes);

%LCCDE Lowpass
LCCDE_Lo = zeros(length(sampleTimes),1);

for n = 2:length(sampleTimes)
    LCCDE_Lo(n) = (1-samplePeriod/tau)*(LCCDE_Lo(n-1)) + samplePeriod/(tau)*inputFunct(n-1);
end

%LCCDE Highpass
LCCDE_Hi = zeros(length(sampleTimes),1);
for n = 2:length(sampleTimes)
    LCCDE_Hi(n) = inputFunct(n) - inputFunct(n-1) - (LCCDE_Hi(n-1) * ((samplePeriod/tau) -1));
end

 % Create subplots
 %plot lowpass
    figure;
    hold on
    subplot(3, 1, 1);
    plot(sampleTimes, abs(inputFunct), 'b', sampleTimes, abs(LCCDE_Lo), 'r');
    title('Magnitude');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 2);
    plot(sampleTimes, real(inputFunct), 'b', sampleTimes, real(LCCDE_Lo), 'r');
    title('Real Part');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 3);
    plot(sampleTimes, imag(inputFunct), 'b', sampleTimes, imag(LCCDE_Lo), 'r');
    title('Imaginary Part');
    legend('Input', 'Output');
    xlabel('Time(s)');
    ylabel('Output');
    sgtitle(['Lowpass Exponential Responses: Frequency = ', num2str(freq), 'Hz']);

    hold off
    % plot highpass
    figure;
    hold on
    subplot(3, 1, 1);
    plot(sampleTimes, abs(inputFunct), 'b', sampleTimes, abs(LCCDE_Hi), 'r');
    title('Magnitude');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 2);
    plot(sampleTimes, real(inputFunct), 'b', sampleTimes, real(LCCDE_Hi), 'r');
    title('Real Part');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 3);
    plot(sampleTimes, imag(inputFunct), 'b', sampleTimes, imag(LCCDE_Hi), 'r');
    title('Imaginary Part');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    sgtitle(['Lowpass Exponential Responses: Frequency = ', num2str(freq), 'Hz']);
end
%% Bode Plot