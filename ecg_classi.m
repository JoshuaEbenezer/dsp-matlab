clear;
clc;
close all;

fnames=dir('train/*.mat');
normal = zeros(7,2,38400);
arrhy = zeros(7,2,38400);

val_normal = zeros(3,2,38400);
val_arrhy = zeros(3,2,38400);

test_normal = zeros(4,2,38400);
test_arrhy = zeros(4,2,38400);
    
Fs = 128;
% number of sampling instances
N = 38400;
% frequency base
f = (-N/2:N/2-1)*Fs/N;
% time base
t = 0:1/Fs:5*60-1/Fs;

%% Read the data and normalize it
for i=1:14
    file=fullfile('train',fnames(i).name);
    s = load(file);
    if (i<=7)
        normal(i,:,:)=s.val;
        % Find energy of the signal and normalize
        energy1 = sum(normal(i,1,:).^2);
        energy2 = sum(normal(i,2,:).^2);
        normal(i,1,:)=normal(i,1,:)/sqrt(energy1);
        normal(i,2,:)=normal(i,2,:)/sqrt(energy2);

    else
        arrhy(i-7,:,:)=s.val;
        % Find energy of the signal and normalize
        energy1 = sum(arrhy(i-7,1,:).^2);
        energy2 = sum(arrhy(i-7,2,:).^2);
        arrhy(i-7,1,:)=arrhy(i-7,1,:)/sqrt(energy1);
        arrhy(i-7,2,:)=arrhy(i-7,2,:)/sqrt(energy2);
    end

end

n = 4;
Wn = 1/(Fs/2);
% Zero-Pole-Gain design
[z,p,k] = butter(n,Wn,'low');
sos = zp2sos(z,p,k);
enorm =zeros(7,1);
earrhy = zeros(7,1);

ncor1 = zeros(7,1);
ncor2 = zeros(7,1);
acor1 = zeros(7,1);
acor2 = zeros(7,1);

for i=1:14
    if (i<=7)
        s1=normal(i,1,:);
        s2=normal(i,2,:);
        
        % filter the signals
        y1 = sosfilt(sos,s1);
        y2 = sosfilt(sos,s2);
        
        % Find the energy of the filtered signals
        enorm(i) = enorm(i)+sum(y1.^2);
        enorm(i) = enorm(i)+sum(y2.^2);
        
        % Find auto-correlations of original signals
        ncor1(i) = sum(s1(1,1,1:end-128).*s1(1,1,129:end));
        ncor2(i) = sum(s2(1,1,1:end-128).*s2(1,1,129:end));
    
    else
        s1=arrhy(i-7,1,:);
        s2=arrhy(i-7,2,:);
        
        y1 = sosfilt(sos,s1);
        y2 = sosfilt(sos,s2);
        % Energy of filtered signals
        earrhy(i-7) = earrhy(i-7)+sum(y1.^2);
        earrhy(i-7) = earrhy(i-7)+sum(y2.^2);
        % Auto-correlation of original signals
        acor1(i-7) = sum(s1(1,1,1:end-128).*s1(1,1,129:end));
        acor2(i-7) = sum(s2(1,1,1:end-128).*s2(1,1,129:end));

    end
end

%% ------------------------------ Validation -------------------------- %%
val_names=dir('val/*.mat');
norm_thresh = zeros(4,1);
arrhy_thresh = zeros(4,1);
confidence = zeros(4,1);

norm_measure = mean(enorm)*(abs(mean(ncor1))+abs(mean(ncor2)))/2;
abnorm_measure = mean(earrhy)*(abs(mean(acor1))+abs(mean(acor2)))/2;

for j=1:4
    weight = 0.6+j*0.1; % Weighing the thresholds for best results
    norm_thresh(j) = norm_measure*weight;
    arrhy_thresh(j) = abnorm_measure*weight;
    
    for i=1:6
        file=fullfile('val',val_names(i).name);
        s = load(file);

        if (i<=3)
            val_normal(i,:,:)=s.val;
            % Find energy of the signal and normalize
            [measure1,measure2]=measure(val_normal,i,sos);
              
            norm_true1 = abs(measure1-norm_thresh(j))<abs(measure1-arrhy_thresh(j));
            norm_true2 = abs(measure2-norm_thresh(j))<abs(measure2-arrhy_thresh(j));
            if ( norm_true1 && norm_true2)
                confidence(j)=confidence(j)+2;
            elseif (~norm_true1 && ~norm_true2)
                confidence(j)=confidence(j)-2;
            end            
        else
            val_arrhy(i-3,:,:)=s.val;
            % Find energy of the signal and normalize
            [measure1,measure2]=measure(val_arrhy,i-3,sos);
            
            abnorm_true1 = abs(measure1-norm_thresh(j))>abs(measure1-arrhy_thresh(j));
            abnorm_true2 = abs(measure2-norm_thresh(j))>abs(measure2-arrhy_thresh(j));

            if (abnorm_true1 && abnorm_true2)
                confidence(j)=confidence(j)+2;
            elseif (~abnorm_true1 && ~abnorm_true2)
                confidence(j)=confidence(j)-2;
            end            
            
        end
    end
end
    
max_chk=fliplr(confidence);
[~,I] = max(max_chk);
weight = 0.6+I*0.1;
fin_norm_thresh = norm_measure*weight;
fin_arrhy_thresh = abnorm_measure*weight;

%% -------------------------------- Testing ---------------------------- %%
test_names=dir('test/*.mat');
tp = 0;
tn = 0;
fp=0;
fn=0;
    
for i=1:4
    file=fullfile('test',test_names(i).name);
    s = load(file);

    if (i<=2)
        test_normal(i,:,:)=s.val;
        % Find energy of the signal and normalize
        
        [measure1,measure2]=measure(test_normal,i,sos);

        norm_true1 = abs(measure1-fin_norm_thresh)<abs(measure1-fin_arrhy_thresh);
        norm_true2 = abs(measure2-fin_norm_thresh)<abs(measure2-fin_arrhy_thresh);
        if (norm_true1)
            tp=tp+1;
        else
            fn=fn+1;
        end
        if (norm_true2)
            tp=tp+1;
        else
            fn=fn+1;
        end         
    else
        test_arrhy(i-2,:,:)=s.val;
        
        [measure1,measure2]=measure(test_arrhy,i-2,sos);

        abnorm_true1 = abs(measure1-fin_norm_thresh)>abs(measure1-fin_arrhy_thresh);
        abnorm_true2 = abs(measure1-fin_norm_thresh)>abs(measure1-fin_arrhy_thresh);

        if (abnorm_true1)
            tn=tn+1;
        else
            fp=fp+1;
        end
        if (abnorm_true2)
            tn=tn+1;
        else
            fp=fp+1;
        end

    end
end
precision=tp/(tp+fn)*100;
recall=tn/(tn+fp)*100;

function [measure1,measure2]=measure(input,i,sos)


        energy1 = sum(input(i,1,:).^2);
        energy2 = sum(input(i,2,:).^2);
        input(i,1,:)=input(i,1,:)/sqrt(energy1);
        input(i,2,:)=input(i,2,:)/sqrt(energy2);

        y1 = sosfilt(sos,input(i,1,:));
        y2 = sosfilt(sos,input(i,2,:));   

        energy1 = sum(y1.^2);
        energy2 = sum(y2.^2);
        test_acor1 = sum(input(i,1,1:end-128).*input(i,1,129:end));
        test_acor2 = sum(input(i,2,1:end-128).*input(i,2,129:end));

        measure1=energy1*abs(test_acor1);
        measure2=energy2*abs(test_acor2);

end
