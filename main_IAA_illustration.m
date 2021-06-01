clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
addpath('E:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;


Sig1 = 1*exp(1i*(1*pi*(18*t.^3)*1)+1i*(2*pi*(10*t))); %300t»òÕß150t
Sig2 = 1*exp(1i*(1*pi*(18*t.^3))*1+1i*(2*pi*(25*t))); %300t»òÕß150t
%Sig2 = 1*exp(1i*(1*pi*(30*t.^3))+1i*(1*pi*(20*t))); %300t»òÕß150t

Sig3 = exp(1i*(1*pi*(120*t -18*t.^3*1)));
%Sig1 = 1*exp(1i*(2*pi*(30*t))); %300t»òÕß150t

Sig =1*Sig1+0*Sig3+1*Sig2;
%Sig =1*Sig1;
%Sig=Sig.*hamming(length(Sig)).';
cccc=30*3;
IF_O(:,1)=cccc*t.^2/2;
IF_O(:,2)=-cccc*t.^2/2+100/2;


%Sig=Sig.*([1:128 128:-1:1]);
num=2;
NS=2;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED

%for snr=-10:2:10
    %    Sig=awgn(SigO,30,'measured');
SigO=Sig;
 %      Sig=awgn(SigO,30,'measured');
%p = randperm(length(Sig));

%p=p(1:44);
%  p=[11:32  65:85  100:111];
% 
% % p=[17:32  85:90  100:106];
 % p=[12:23  66:75  100:109];
% % 
p=[12:30 40:55 66:75  100:109];

p=[10:30 41:54 66:80  95:109];

% p=[20:35 60:75 95:115];

Sig(p)=0;
iii=find(Sig~=0);
% ORIGINAL
delta=5;
alpha = 5;
mean(abs(SigO-Sig))

[Sig_out1,findex] = FAST_IF_Recover(Sig,length(Sig)/(2)-1, num, 2,100,0,0,iii);
[Sig_out,findex] = FAST_IF_Recover_iterative(Sig,length(Sig)/(2)-1, num, 2,100,0,0,iii);

%[Sig_out,findex] = FAST_IF_Recover1(Sig,length(Sig)/(2)-1, num, 3,100,0,0,iii);

figure; plot(real(SigO),'r:','Linewidth',2);
 hold on; plot(real(Sig),'k','Linewidth',2);
%plot(real(SigO))
hold on;plot(real(Sig_out1),'b--','Linewidth',2);
hold on;plot(real(Sig_out),'g-.','Linewidth',2);

xlabel('Samples')
ylabel('Amplitude')
legend('The Original Signal','Sparsely Sampled Signal','The Reconstructed Signal (Single Iteration)','The Reconstructed Signal (3 Iterations)');

 I1=HTFD_new1(SigO,3,8,64);
 I2=HTFD_new1(Sig,3,8,64);
 I3=HTFD_new1(Sig_out1,3,8,64);
 I4=HTFD_new1(Sig_out,3,8,64);

figure; 
imagesc(I1)
figure;
imagesc(I2)
figure;
imagesc(I3)
figure;
imagesc(I4);
