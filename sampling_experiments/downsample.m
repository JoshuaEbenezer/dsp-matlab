%% Script file for Exp 1 in DSP lab
clc;
clear;
close all;

%% --------------------------------- Part D ---------------------------- %%
% sampling freq (this was changed for part B)
Fs=12e3;
% time interval
Ts = 1/Fs;
% original time interval
Tsorig = 1/96e3;
% original frequency of sampling
Fsorig = 1/Tsorig;
% length of sampling base
L = 500;
% time base original
torig = (0:(Fsorig/Fs)*L)*Tsorig;

% sampling time base
t = (0:L)*Ts;

% generate signal
xorig = 10*cos(2*pi*1e3*torig)+6*cos(2*pi*2e3*torig)+2*cos(2*pi*4*1e3*torig);

figure;
plot(torig(1:100),xorig(1:100));
xlabel('time(s)');
ylabel('Original signal');

% sampled signal
x = 10*cos(2*pi*1e3*t)+6*cos(2*pi*2e3*t)+2*cos(2*pi*4*1e3*t);

% frequency base for sampled signal
f1 = (-L/2:L/2)*Fs/L;

% generate fft
fft_x = fftshift(fft(x));
%find amplitudes
amp_fft_x1 = abs(fft_x/L);

figure;
stem(f1,amp_fft_x1);
title('Amplitude spectrum of sampled signal');
xlabel('f(Hz)');
ylabel('Amplitude');

% downsample the signal by a factor of 5
M = 5;
down_x = downsample(x,M);
t_down = (0:L/M)*(M*Ts);

%% Perform sample and hold on the sampled signal

% zeroth order interpolation
samp_hold = interp1(t_down,down_x,torig,'previous');
figure;
plot(torig,samp_hold);
title('Zero order hold sampling without anti-aliasing');
xlabel('time(s)');
ylabel('Amplitude');

% cutoff frequency
fc = 6e3/(Fsorig/2);
% LPF
[z,p,k] = butter(4,fc,'low');
sos = zp2sos(z,p,k);
% recover signal via lpf
zero_order_recover = sosfilt(sos,samp_hold);

figure;
plot(torig,zero_order_recover);
hold on;
plot(torig,xorig,':r','LineWidth',2);
title('Zero order hold recovery without anti-aliasing');
xlabel('time(s)');
ylabel('Amplitude');

% length of the downsampled signal
L_down = L/M;
% frequency of the downsampled signal
Fs_down = Fs/M;

% generate fft
fft_x1 = fftshift(fft(down_x));
%find amplitudes
amp_fft_x1 = abs(fft_x1/L_down);

% frequency base for the downsampled signal
f_down = (-L_down/2:L_down/2)*(Fs_down/L_down);

figure;
stem(f_down,amp_fft_x1);
title('Amplitude spectrum of downsampled signal without anti-aliasing');
xlabel('f(Hz)');
ylabel('Amplitude');

%% Anti aliasing filter

% cutoff frequency
fc_antialias = 1/M;
% LPF
[z,p,k] = butter(4,fc_antialias,'low');
sos = zp2sos(z,p,k);
% apply anti aliasing filter
antialias = sosfilt(sos,x);

% downsample output of anti aliasing filter
down_anti = downsample(antialias,M);
 
%% FFT of downsampled output of anti aliasing filter

% frequency base for the downsampled signal
f_down = (-L_down/2:L_down/2)*(Fs_down/L_down);

% generate fft
fft_x1 = fftshift(fft(down_anti));
%find amplitudes
amp_fft_x1 = abs(fft_x1/L_down);

figure;
stem(f_down,amp_fft_x1);
title('Amplitude spectrum of downsampled signal after anti-aliasing');
xlabel('f(Hz)');
ylabel('Amplitude');


% zeroth order interpolation
anti_alias_reco = interp1(t_down,down_anti,torig,'previous');

% generate fft
fft_x1 = fftshift(fft(anti_alias_reco));
%find amplitudes
amp_fft_x1 = abs(fft_x1/length(anti_alias_reco));

figure;
plot(torig,anti_alias_reco);
title('Zero order hold recovery after anti-aliasing');
xlabel('time(s)');
ylabel('Amplitude');

figure;
plot(torig,zero_order_recover);
hold on;
plot(torig,xorig,':r','LineWidth',2);
title('Signal recovered from sample and hold after downsampling and anti-aliasing');
xlabel('time(s)');
ylabel('Amplitude');

set(findall(gcf,'-property','FontSize'),'FontSize',24);
