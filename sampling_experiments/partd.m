%% ------------------Script file: Part d) of EXP 1 (DSP lab) ----------%%

clc;
clear;
close all;

% sampling freq 
Fs=12e3;
% time interval
Ts = 1/Fs;
% original analog time interval
Tsorig = 1/48e3;
% 'analog' sampling freq
Fsorig = 1/Tsorig;
% length of sampling base
L = 1000;
% time base original
torig = (0:(Fsorig/Fs)*L-1)*Tsorig;

% sampling time base
t = (0:L-1)*Ts;

% generate original signal
xorig = 10*cos(2*pi*1e3*torig)+6*cos(2*pi*6e3*torig)+2*cos(2*pi*4*1e3*torig);
% sampled signal
x = 10*cos(2*pi*1e3*t)+6*cos(2*pi*6e3*t)+2*cos(2*pi*4*1e3*t);

%% --------------------------FFT of sampled signal --------------------- %%
% frequency base for sampled signal
f1 = (-L/2:L/2-1)*Fs/L;

% generate fft
fft_x = fftshift(fft(x));
%find amplitudes
amp_fft_x1 = abs(fft_x/L);

figure;
stem(f1,amp_fft_x1);
title('Amplitude spectrum of sampled signal');
xlabel('f(Hz)');
ylabel('Amplitude');

%% ------------------------ upsample the signal ------------------------ %%
x_upsampl = upsample(x,2);

%% ----------------------- Display signal FFT before LPF ----------------%%
% length of upsampled signal
L_up = L*2;
% frequency of sampling
F_up = Fs*2;
% frequency base
f_up = (-L_up/2:L_up/2-1)*(F_up/L_up);

% generate fft
fft_x = fftshift(fft(x_upsampl));
%find amplitudes
amp_fft_x1 = abs(fft_x/L_up);

figure;
stem(f_up,amp_fft_x1);
title('Amplitude spectrum of upsampled signal before LPF');
xlabel('f(Hz)');
ylabel('Amplitude');

%% ------------------------ Low pass filter ---------------------------- %%

% find filter coefficients for 4th order LPF
lpf = fir1(4,6e3/F_up,'low'); % FIR clock is 24 kHz

% pass signal through LPF
out = filter(lpf,1,x_upsampl);

% generate fft
fft_x = fftshift(fft(out));
%find amplitudes
amp_fft_x1 = abs(fft_x/L_up);

figure;
stem(f_up,amp_fft_x1);
title('Amplitude spectrum of upsampled signal after LPF');
xlabel('f(Hz)');
ylabel('Amplitude');

%% ------------------------------ Comparison --------------------------- %%

% sample original signal at 24 kHz
t24sampl = (1:24e3)/24e3;

x24sampl = 10*cos(2*pi*1e3*t24sampl)+6*cos(2*pi*6e3*t24sampl)+2*cos(2*pi*4*1e3*t24sampl);
 
% plot original signal sampled at 24 kHz and output of LPF
figure;
plot(t24sampl(1:100),x24sampl(1:100));
hold on;
plot(t24sampl(1:100),out(1:100),'r');
xlabel('time(s)');
ylabel('Amplitude');
title('Output of LPF and original signal');
