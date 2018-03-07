clc;
clear;
close all;

% sampling freq 1
Fs1=20e3;
% time interval for sampling
Ts1 = 1/Fs1;
% time base for sampling
tsamp1 = (0:1000)*Ts1;
% last time instance
tmax = 1000*Ts1;

% original time base
t = 0:1/(100e3):tmax;
% generate square wave
xsq_orig = square(2*pi*1e3*t); 
% sampled square wave
xsq1 = square(2*pi*1e3*tsamp1);
% no of points for FFT
N = 256;
% frequency base
f1 = (-N/2:N/2-1)*Fs1/N;
% find fft
fft_xsq1 = fftshift(fft(xsq1,N));
% find magnitudes of fft
abs_fft1 = abs(fft_xsq1/N);

% plot original signal
figure;
plot(t(1:500),xsq_orig(1:500));
title('Original signal (truncated)');
xlabel('time(s)');
ylabel('Original signal');

figure;
stem(tsamp1(1:100),xsq1(1:100));
title('Sampled signal (truncated)');
xlabel('time(s)');
ylabel('Sampled signal');

% plot
figure;
plot(f1,abs_fft1);
title('Amplitude spectrum of square wave for N=256, Fs=20 kHz');
xlabel('f(Hz)');
ylabel('Amplitude');

fs = zeros(1,50);
ts = zeros(1,50);
mean_sq_er = zeros(1,50);
for i=1:50
    fs(i) = i*1000; %k*fo
    ts(i) = 1/fs(i);
    tsamp = 0:ts(i):tmax;
    % sampled square wave
    xsq = square(2*pi*1e3*tsamp);
    
    % reconstruct from the sampled wave
    recons = interp1(tsamp,xsq,t);
    % find mean square error
    mean_sq_er(i) = sum((recons-xsq_orig).^2)/length(t);
end
figure;
plot(fs/1000,mean_sq_er);
title('Mean square error vs sampling frequency');
xlabel('Sampling frequency (kHz)');
ylabel('Mean Square error');






