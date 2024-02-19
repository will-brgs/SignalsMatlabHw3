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

sampleTimes = 0:samplePeriod:15*(1/freq);
%inputFunct = zeros(1,length(timePoints));
inputFunct = exp(1j*angularFreq*sampleTimes);

%LCCDE Lowpass

lsim_Lo = lsim(1/tau , [1 ,1/tau], inputFunct,sampleTimes);

%LCCDE Highpass

lsim_Hi = lsim([1 0] , [1 ,1/tau], inputFunct,sampleTimes);

 % Create subplots
 %plot lowpass
    figure;
    hold on
    subplot(3, 1, 1);
    plot(sampleTimes, abs(inputFunct), 'b', sampleTimes, abs(lsim_Lo), 'r');
    title('Magnitude');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 2);
    plot(sampleTimes, real(inputFunct), 'b', sampleTimes, real(lsim_Lo), 'r');
    title('Real Part');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 3);
    plot(sampleTimes, imag(inputFunct), 'b', sampleTimes, imag(lsim_Lo), 'r');
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
    plot(sampleTimes, abs(inputFunct), 'b', sampleTimes, abs(lsim_Hi), 'r');
    title('Magnitude');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 2);
    plot(sampleTimes, real(inputFunct), 'b', sampleTimes, real(lsim_Hi), 'r');
    title('Real Part');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    
    subplot(3, 1, 3);
    plot(sampleTimes, imag(inputFunct), 'b', sampleTimes, imag(lsim_Hi), 'r');
    title('Imaginary Part');
    xlabel('Time(s)');
    ylabel('Output');
    legend('Input', 'Output');
    sgtitle(['Highpass Exponential Responses: Frequency = ', num2str(freq), 'Hz']);
end
%% Part 2 Bode Plot
freqRange = logspace(1,4,100);
sampleTimes = 0:samplePeriod:2; %Sim to 2s for bettter steady state

H_Hi = zeros(length(freqRange),1);
H_Lo = zeros(length(freqRange),1);

for i = 1:length(freqRange)
    angularFreq = freqRange(i);
    inputFunct = exp(1j*angularFreq*sampleTimes);

    %Lsim Lowpass
    lsim_Lo = lsim(1/tau , [1 ,1/tau], inputFunct,sampleTimes);

    %Lsim Highpass
    lsim_Hi = lsim([1 0] , [1 ,1/tau], inputFunct,sampleTimes);

    % Compute H, assume steady state at index 660
    H_Lo(i) =  lsim_Lo(end) / inputFunct(end);
    H_Hi(i) =  lsim_Hi(end) / inputFunct(end);
    % Compute magnitude in dB and phase normalized by pi

end

%Calc for Highpass
mag_Hi = 20*log10(abs(H_Hi));
phase_Hi = angle((H_Hi)/pi);

%Calc for Lowpass
mag_Lo = 20*log10(abs(H_Lo));
phase_Lo = angle((H_Lo)/pi);

%Plot
figure;
hold on
subplot(4, 1, 1);
semilogx(freqRange, mag_Lo, LineWidth=1.5);
title('Lowpass Magnitude');
xlabel('Frequency (Hz)');
ylabel('Output (dB)');
    
subplot(4, 1, 2);
semilogx(freqRange, mag_Hi, LineWidth=1.5);
title('HighPass Magnitude');
xlabel('Frequency (Hz)');
ylabel('Output (dB)');
    
subplot(4, 1, 3);
semilogx(freqRange, phase_Lo, LineWidth=1.5);
title('LowPass Phase');
xlabel('Frequency (Hz)');
ylabel('Output Phase (Radians)');

subplot(4, 1, 4);
semilogx(freqRange, phase_Hi, LineWidth=1.5);
title('HighPass Phase');
xlabel('Frequency (Hz)');
ylabel('Output Phase (Radians)');
sgtitle('Bode Plot Outputs of High and Lowpass Filters with Varying Frequency');
hold off

%% Part 3: DC to AC Converter
% Define parameters
sampleFreq = 10e3;
samplePeriod = 1/sampleFreq;
frequency = 60;      % Desired frequency in Hz

cutoff = 160;
C = 10e-6;
R = 1/(2*pi*cutoff*C);
tau = R*C;

sampleTimes = 0:samplePeriod:20 *(1/frequency);

% Generate square wave
sq = square(2*pi*frequency*sampleTimes);
inputFunct = sq;

%Filter Coeffs
a = [1, samplePeriod/tau-1];
b = samplePeriod/tau;

% Pass signal through RC Lowpass filter 6 separate times
output = zeros(size(inputFunct));
for i = 1:6
output = filter(b, a, inputFunct);

inputFunct = output;
end
 
% Plot the output
figure()
hold on
plot(sampleTimes, sq);
plot(sampleTimes, output);
xlabel('Time (s)')
ylabel('Output (V)')
title('DC to AC Converter')
legend('Input', 'Output');
hold off
%% Part 3C: Power Efficiency
Pin = (1/samplePeriod) * sum(sq(end-sampleFreq*(1/frequency):end).^2);
% 218 ish datapoints per period

Pout = (1/samplePeriod) * sum(output(end-sampleFreq*(1/frequency):end).^2);

efficiency = Pout/Pin;
fprintf('Efficiency: %f\n', efficiency)

%% Part 3B : Bode Plot
freqRange = logspace(1,4,100);

H = zeros(length(freqRange),1);
for i = 1:length(freqRange)
    frequency = freqRange(i);
    sq = square(2*pi*frequency*sampleTimes);
    inputFunct = sq;

    % Pass signal through RC Lowpass filter 3 separate times
    output = zeros(size(inputFunct));
    for j = 1:6
    output = filter(b, a, inputFunct);
    
    inputFunct = output;
    end

    % Compute H, assume steady state at end
    H(i) =  output(end) / sq(end);
    % Compute magnitude in dB and phase normalized by pi

end

%Calc for Lowpass
mag = 20*log10(abs(H));
phase = angle((H)/pi);

%Plot
figure;
hold on
subplot(2, 1, 1);
semilogx(freqRange, mag, LineWidth=1.5);
title('Bode Magnitude');
xlabel('Frequency (Hz)');
ylabel('Output (dB)');
    
subplot(2, 1, 2);
semilogx(freqRange, phase, LineWidth=1.5);
title('Bode Phase');
xlabel('Frequency (Hz)');
ylabel('Output');
sgtitle('Bode Plot Output of DC to AC Converter with Varying Frequency');
hold off
