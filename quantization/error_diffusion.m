%% Script for error diffusion quantization
clc;
clear;
close all;

img = imread('vict.jpg');
gray_img = rgb2gray(img);
quant_img = zeros(size(gray_img));

[r,c]=size(gray_img);
error = double(gray_img);

for i=1:r
    if (mod(i,2)==1)    % if in odd row then go from 1st col to end
        for j=1:c       % iterate through columns
            if(error(i,j)<=127) % threshold - half grey
                quant_img(i,j)=0; 
                if (j~=c && i~=r) % if not at the end of row/col
                    error(i,j+1)=error(i,j+1)+error(i,j); 
                elseif(i~=r) % if not at end of row
                    error(i+1,j)=error(i+1,j)+error(i,j);
                end % otherwise no diffusion required
            else
                quant_img(i,j)=255;
                if (j~=c && i~=r)
                    error(i,j+1)=error(i,j+1)+error(i,j)-255; % diffuse error
                elseif(i~=r)
                    error(i+1,j)=error(i+1,j)+error(i,j)-255;
                end
            end
        end
    else      % if in even row go from last column to first
        for j=c:-1:1
            if(error(i,j)<=127)
                quant_img(i,j)=0;
                if (j~=1 && i~=r)
                    error(i,j-1)=error(i,j-1)+error(i,j);
                elseif(i~=r)
                    error(i+1,j)=error(i+1,j)+error(i,j);
                end
            else
                quant_img(i,j)=255;
                if (j~=1 && i~=r)
                    error(i,j-1)=error(i,j-1)+error(i,j)-255;
                elseif(i~=r)
                    error(i+1,j)=error(i+1,j)+error(i,j)-255;
                end
            end
        end
    end
end

