%% Script to study quantization
clc;
close all;
clear;

img = imread('vict.jpg');
gray_img = rgb2gray(img);
img_hist = imhist(gray_img);
mse = zeros(1,7);
%store_b = zeros(10,4);
for bits=2 % number of bits = log(n)
    n = 2^bits; % n - number of bins
    b = zeros(1,n);

    for i=1:n
        b(i) =255*i/n;
    end
    itermax = 10;   % ten iterations
    error = zeros(itermax,n);
    mse_loop = zeros(itermax,1);
    for iter = 1:itermax
    exp_val = zeros(1,n);   % to store expected value E(x|Ri)
    denom = zeros(1,n);
       for pixel=1:256  %iterate through pixels
           for i=1:n    % find which region it belongs to
               if(pixel-1<=b(i)) % defining each region Ri (bin)
                   exp_val(i)=exp_val(i)+(pixel-1)*img_hist(pixel); %finding expected value E(x|Ri)
                   denom(i) = denom(i)+img_hist(pixel); % p(Ri) being found to divide at the end
                   break;
               end
           end
       end
       exp_XR=exp_val./denom;   % expression for E(x|R)
       %store_b(iter,:) = b; 
       new_b = [(exp_XR(1:end-1)+exp_XR(2:end))/2, 255];
       error(iter,:) = new_b-b;
       mse_loop(iter) = mean(error(iter,:).^2);
       b=new_b; % update b

    end
    % quantize according to b
    quant_img = zeros(size(gray_img),'uint8');
    quant_img(gray_img<=b(1))=uint8(b(1)/2);
    for i=1:n-1
        quant_img(gray_img<=b(i+1) & gray_img>b(i))=uint8((b(i)+b(i+1))/2);
    end
    mse(bits-1) = mean((quant_img(:)-gray_img(:)).^2);
end
% figure;
% plot(1:itermax,mse_loop);
% xlabel('Number of iterations','fontsize',20);
% ylabel('Mean square error','fontsize',20);
% title('MSE vs number of iterations','fontsize',24);
figure;
plot(2:8,mse);
xlabel('Number of bits','fontsize',20);
ylabel('Mean square error','fontsize',20);
title('MSE vs number of bits','fontsize',24);