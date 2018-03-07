%% Script to study quantization
clc;
close all;
clear;
% Sampling frequency
Fs = 1000;
% Sampling time interval
Ts = 1/Fs;
% time base
t = (0:4000-1)*Ts;
% sine wave
x = sin(2*pi*5*t);

%% Uniform quantization
delta = 1/2^11;
partition = -2048*delta+delta/2:delta:2048*delta-delta/2;
codebook=-2048*delta:delta:2048*delta;

[~,unif_quant] = quantiz(x,partition,codebook);
precision = 12;
unif_binary = zeros(4000,12);
%% Convert to binary
for m = 1:4000
    val = abs(unif_quant(1,m));
    for i=1:precision-1
        unif_binary(m,i+1) = floor(val*2);
        val = val*2-floor(val*2);
    end
    if (sign(unif_quant(m))==-1)
        unif_binary(m,1)=1;
    end
end

%% truncate to 8 bits
truncated = unif_binary(:,1:8);
dec_trunc = zeros(4000,1);
error = zeros(4000,1);
places = pow2(-(1:7));
for m=1:4000
    absval = truncated(m,2:end);
    dec_trunc(m) = (-1)^truncated(m,1) * sum(truncated(m,2:end).*places,2);
    error(m) = sum(unif_quant(m)-dec_trunc(m));
end
%% FFT
N = 4000;
% frequency base
f = (-N/2:N/2-1)*Fs/N;
% generate fft
fft_xu = fftshift(fft(dec_trunc,N));
%find amplitudes
amp_fft_xu = abs(fft_xu/N);

figure;
stem(f,amp_fft_xu);
title('Amplitude spectrum for truncated signal','fontsize',24);
xlabel('f(Hz)','fontsize',24);
ylabel('Amplitude','fontsize',24);

mse = mean(error.^2);

