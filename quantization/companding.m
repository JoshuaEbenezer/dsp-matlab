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

u = [7,15,31,127,255];
A = [1,2,5,10,50,87.6];
u_compand=zeros(5,4000);
a_compand=zeros(6,4000);
val = zeros(1,4000);

%% plot u and A law transfer functions
dummy = linspace(-1,1,4000);
utrans = sign(dummy).*log(1+u(3).*abs(dummy))./log(1+u(3));

k1 = find(abs(dummy)<1/A(6));
val(k1) = A(6)*abs(dummy(k1))/(1+log(A(6)));

k2 = find(abs(dummy)>1/A(6));
val(k2)= (1+log(A(6)*abs(dummy(k2))))/(1+log(A(6)));

atrans = sign(dummy).*val;

figure;
plot(dummy,atrans);
hold on;
plot(dummy,utrans,'r');
ylabel('Output','fontsize',20);
xlabel('Input','fontsize',20);
legend('A-law','u-law');
title(['Transfer functions for A = ' num2str(A(6)) ' and u = ' num2str(u(3))],'fontsize',24 );


%% u-law and A law transforms
for index=1:5
    u_compand(index,:) = sign(x).*log(1+u(index).*abs(x))./log(1+u(index));
    
    k1 = find(abs(x)<1/A(index));
    val(k1) = A(index)*abs(x(k1))/(1+log(A(index)));
   
    k2 = find(abs(x)>1/A(index));
    val(k2)= (1+log(A(index)*abs(x(k2))))/(1+log(A(index)));
    
    a_compand(index,:) = sign(x).*val;
end
index=index+1;
k1 = find(abs(x)<1/A(index));
val(k1) = A(index)*abs(x(k1))/(1+log(A(index)));

k2 = find(abs(x)>1/A(index));
val(k2)= (1+log(A(index)*abs(x(k2))))/(1+log(A(index)));
a_compand(6,:) = sign(x).*val;


%% quantization
delta = 1/2^7;
partition = -128*delta+delta/2:delta:127*delta-delta/2;
codebook=-128*delta:delta:127*delta;
a_par_index = zeros(6,4000);
a_quant = zeros(6,4000);

u_par_index = zeros(6,4000);
u_quant = zeros(5,4000);

aerror = zeros(6,4000);
uerror = zeros(5,4000);

for i=1:5
    [a_par_index(i,:),a_quant(i,:)] = quantiz(a_compand(i,:),partition,codebook);
    [u_par_index(i,:),u_quant(i,:)] = quantiz(u_compand(i,:),partition,codebook);
    
    aerror(i,:) = a_quant(i,:)-x;
    uerror(i,:) = u_quant(i,:)-x;
%     figure;
%     stem(a_quant(i,1:200));
%     hold on;
%     stem(u_quant(i,1:200),'r');
%     stem(x(1:200),'k');
%     ylabel('Amplitude/Quantized value','fontsize',20);
%     xlabel('Samples','fontsize',20);
%     legend('A-law','u-law','Original signal');
%     title(['A = ' num2str(A(i)) ' u = ' num2str(u(i))],'fontsize',24 );

%     %% FFT
%     N = 4000;
%     % frequency base
%     f = (-N/2:N/2-1)*Fs/N;
%     % generate fft
%     fft_xu = fftshift(fft(u_quant(i,:),N));
%     %find amplitudes
%     amp_fft_xu = abs(fft_xu/N);
% 
%     figure;
%     stem(f,amp_fft_xu);
%     title(['Amplitude spectrum for u = ', num2str(u(i))]);
%     xlabel('f(Hz)');
%     ylabel('Amplitude');

%     % generate fft
%     fft_xa = fftshift(fft(a_quant(i,:),N));
%     %find amplitudes
%     amp_fft_xa = abs(fft_xa/N);
% 
%     figure;
%     stem(f,amp_fft_xa);
%     title(['Amplitude spectrum for A = ', num2str(A(i))]);
%     xlabel('f(Hz)');
%     ylabel('Amplitude');
end
[a_par_index(6,:),a_quant(6,:)] = quantiz(a_compand(6,:),partition,codebook);

aerror(6,:) = a_quant(6,:)-x;

% figure;
% stem(a_quant(6,1:200));
% hold on;
% stem(x(1:200),'k');
% ylabel('Amplitude/Quantized value','fontsize',20);
% xlabel('Samples');
% legend('A-law','Original signal');
% title(['A = ' num2str(A(6))],'fontsize',24);


% % generate fft
% fft_xa = fftshift(fft(a_quant(6,:),N));
% %find amplitudes
% amp_fft_xa = abs(fft_xa/N);
% 
% figure;
% stem(f,amp_fft_xa);
% title(['Amplitude spectrum for A = ' num2str(A(6))]);
% xlabel('f(Hz)');
% ylabel('Amplitude');

%% Uniform quantization
[~,unif_quant] = quantiz(x,partition,codebook);
uniferror = unif_quant-x;
% find mean errors
mean_aerror = mean(abs(aerror),2);
mean_uerror = mean(abs(uerror),2);
mean_uniferror = mean(abs(uniferror),2);

figure;
stem(unif_quant(1:200));
hold on;
stem(x(1:200),'k');
ylabel('Amplitude/Quantized value','fontsize',20);
xlabel('Samples','fontsize',20);
legend('Uniform','Original signal');
title('Uniform quantization','fontsize',24);