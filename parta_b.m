%% Script file for Exp 1 in DSP lab
clc;
clear;
close all;

%% ----------------------------- Part A and B -------------------------- %%
% sampling freq (this was changed for part B)
Fs=12e3;
% time interval
Ts = 1/Fs;
% original time interval
Tsorig = 1/96e3;
Fsorig = 1/Tsorig;
% time base original
torig = (0:4000)*Tsorig;
% sampling time base
t = (0:500)*Ts;

% generate signal
xorig = 10*cos(2*pi*1e3*torig)+6*cos(2*pi*2e3*torig)+2*cos(2*pi*4*1e3*torig);

figure;
plot(torig(1:100),xorig(1:100));
xlabel('time(s)');
ylabel('Original signal');

% sampled signal
x = 10*cos(2*pi*1e3*t)+6*cos(2*pi*2e3*t)+2*cos(2*pi*4*1e3*t);

sn = round(50/(Fs*Tsorig)); % x axis limit for sampled display

figure;
stem(t(1:sn),x(1:sn));
xlabel('time(s)');
ylabel('Sampled signal');

%% Case 1 N = 64
% number of sampling instances
N1 = 64;
% frequency base
f1 = (-N1/2:N1/2-1)*Fs/N1;

% generate fft
fft_x1 = fftshift(fft(x,N1));
%find amplitudes
amp_fft_x1 = abs(fft_x1/N1);

figure;
stem(f1,amp_fft_x1);
title('Amplitude spectrum for N=64');
xlabel('f(Hz)');
ylabel('Amplitude');

%% Case 2 N = 128
% number of sampling instances
N2 = 128;

% frequency base
f2 = (-N2/2:N2/2-1)*Fs/N2;

% generate fft
fft_x2 = fftshift(fft(x,N2));
%find amplitudes
amp_fft_x2 = abs(fft_x2/N2);

figure;
stem(f2,amp_fft_x2);
title('Amplitude spectrum for N=128');
xlabel('f(Hz)');
ylabel('Amplitude');

%% Case 3 N = 256
% number of sampling instances
N3 = 256;
% frequency base
f3 = (-N3/2:N3/2-1)*Fs/N3;
% generate fft
fft_x3 = fftshift(fft(x,N3));
%find amplitudes
amp_fft_x3 = abs(fft_x3/N3);

figure;
stem(f3,amp_fft_x3);
title('Amplitude spectrum for N=256');
xlabel('f(Hz)');
ylabel('Amplitude');

%% Perform sample and hold on the sampled signal

% zeroth order interpolation
samp_hold = interp1(t,x,torig,'previous');
figure;
plot(torig,samp_hold);
title('Zero order hold sampling');
xlabel('time(s)');
ylabel('Sample and hold');

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
title('Signal recovered from sample and hold');
xlabel('time(s)');
ylabel('Recovered signal');

%% First order interpolation

% first order interpolation
lin_inter = interp1(t,x,torig);
figure;
plot(t,x,'o',torig,lin_inter,':.','LineWidth',2);
title('First order hold sampling');
xlabel('time(s)');
ylabel('Amplitude');

figure;
plot(torig,lin_inter);
hold on;
plot(torig,xorig,':r','LineWidth',2);
title('Signal recovered from linear interpolation');
xlabel('time(s)');
ylabel('Recovered signal');
