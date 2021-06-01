clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;

Sig=bat_signal.';

%load seizure;
Sig=hilbert(Sig);

%IF_O(:,3)=cccc*t.^2/2+20/2;
%IF_O(:,4)=-cccc*t.^2/2+115/2;

SigO=Sig;

%IF_O(:,3)=90*t.^2/2+15;
WN=64;
wind_step=32;

%Sig=Sig.*([1:128 128:-1:1]);
num=3;
NS=100;
% HADTFD BASED
iiii=0;
%for snr=-10:2:10
% for N_S=5:5:15
      for N_S=10:10:40

         iiii=iiii+1;
    for k1=1:NS
        Sig=SigO;
        p=[];
        for i=2:7
        pp = 50*(i-1)+ randperm(50-N_S-1,1);
        p1=pp:1:pp+N_S;
        p=[ p p1];
        end
        Sig(p)=0;
        [NA]=find(Sig~=0);
        for kkkkk=0:3
            
            % ORIGINAL
            delta=5;
            alpha = 5;
            if kkkkk==0   %ADTFD+VITERBI
                
                             [ext_sig,findex] = sparse_reconstruction_FAST_IF(Sig, num,11,0.1,p,5,121);
                
            elseif kkkkk==1 %the new algorithm
                
             %   [ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,121, num, 2,100,0,0,NA,7);
                  [ext_sig_iaa,findex] = FAST_IF_Recover_bat(Sig,121, num+2, 10,100,0.1,0,NA);
              %  [ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,121, num, 2,50,0,0,NA,5);
            elseif  kkkkk==2
               [ext_sig] = STFT_RECONSTRUCTION(Sig,WN,wind_step,NA);
            else 
                            y=GradRec(Sig,p);

            end
            
            if kkkkk==0
                                mse_FAST_IF_TF_FILTER(k1)=mean(abs(ext_sig-SigO));

            elseif kkkkk==1
                                mse_FAST_IF_IAA(k1)=mean(abs(ext_sig_iaa-SigO));

            elseif kkkkk==2
                mse_STFT(k1)=mean(abs(ext_sig-SigO));
            else
                mse_GD(k1)=mean(abs(y-SigO));

            end
            
            
            
        end
        
        
    end
    mse_IAA(iiii)=mean(mse_FAST_IF_IAA)
    mse_TF(iiii)=mean(mse_FAST_IF_TF_FILTER)
    mse_ST(iiii)=mean(mse_STFT)
        mse_GDD(iiii)=mean(mse_GD)

 end
 %plot(real(Sig))
 %hold on; plot(real(SigO),'r:');
 %hold on; plot(real(ext_sig_iccd),'k-');
 figure;

plot(10:10:40,mse_IAA(1:4),'k','linewidth',3);
hold on;
plot(10:10:40,mse_ST(1:4),'r','linewidth',3);
hold on;
plot(10:10:40,mse_TF(1:4),'g:','linewidth',3);
hold on;
plot(10:10:40,mse_GDD(1:4),'b','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
legend('The Proposed Method','Modified OMP','TF filtering','Gradient Descent');
set(gca,'fontsize', 36)
%  
%  mse_IAA =
% 
%     0.0034    0.0030    0.0021    0.0015
% 
% 
% mse_TF =
% 
%     0.0046    0.0079    0.0146    0.0261
% 
% 
% mse_ST =
% 
%     0.0088    0.0176    0.0301    0.0403

