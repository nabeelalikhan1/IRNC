

function [sig_out,fidexmult] = FAST_IF_Recover_bat(Sig,win_length, num, delta,L,thr,Thr,iiii)
% SET OF FRACTIONAL WINDOWS
% Input Parameters:
% Sig: Analytical Signal computed using Hilbert Transform
% wind_length: Length of the analysis window
% num: Maximum number of iterations that can be set equal to the length of number of components
% delta: Delta F , i.e maximum possible deviation between the two consecutive sample of the signal.
% L: 2*L+a are the number of quantization levels for the computation of fractional Fourier windows
% thr: Maximum possible ratio between the energy of the strongest component and the weakest component. It is used to stop the algorithm when IFs of all components have been estimated.
% Thr: Maximum possible ratio between the amplitude of the strongest TF point of a given component and the weakest TF point. It is used to find boundary of the IF curve
%Output:
%fidexmult: IFs of the signal component
%TFC: 3D representation of IF curve in time-frequency-chirp domain
w=gausswin(win_length,1);
l=0:win_length-1;
sig_out=zeros(size(Sig));
Siggg=Sig;
for iii=1:length(Sig)
    WW(iii,:)=exp(-1i*(iii)*2*pi*l/length(Sig));
end


indexs=(find(Sig==0));


%A1=AA(:,iiii).';
%A2=AA(:,indexs).';


e_max=0;
%tic;
i=0;
window_rot=zeros(2*L+1,win_length);
for k=-L+1:1:L-1
    i=i+1;
    window_rot(i,:)=frft(w,0.85* k/L);%fracft(w,0.95* k/L);%0.05
end
%save('window_rot','window_rot');
%load('window_rot');


v=zeros(1,2*L+1);
index=v;

for iii=1:num
    Sig_extended=[zeros(1,floor(win_length/2)) Sig zeros(1,floor(win_length/2))];
    
    Siga=filter(ones(1,win_length),1,abs(Sig));
    [~,t_start]=max(Siga(floor(win_length/2)+1:end-floor(win_length/2)));
    t_start=t_start(1)+floor(win_length/2);
    
    for i=1:2*L+1
        FF=abs(fft(Sig(t_start-floor(win_length/2):t_start+floor(win_length/2)).*window_rot(i,:),length(Sig)));
        
        [v(i),index(i)]=max(FF(1:end/2));
    end
    [v_m,ind]=max(v);
    
    freq_start=index(ind);
    frac_wind_index_start=ind;
    v_oldd=v_m;
    
    
    IF=zeros(1,length(Sig))-1;
    IF(t_start)=freq_start-1;
    
    
    clear v;
    
    
    for it=1:2
        
        f0=freq_start;
        frac_wind_index=frac_wind_index_start;
        t_index=t_start;
        while and(t_index>1,t_index<length(Sig))
            
            if it==1
                t_index=t_index+1;
            else
                t_index=t_index-1;
            end
          
            k=f0-delta:1:f0+delta;
            k(k>length(Sig)/2-1)=length(Sig)/2-1;
            k(k<=0)=1;
            if frac_wind_index<2
                frac_wind_index=2;
            elseif frac_wind_index>2*L-1
                
                frac_wind_index=2*L-1;
            end
            %   for i=frac_wind_index-1:frac_wind_index+1    % FOR ALL WINDOWS
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(frac_wind_index-1,:)); %WINDOWED SIGNAL
            V=abs(WW(k,:)*w_signal.');
            [v(1),indexx(1)]=max(V);
            
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(frac_wind_index,:)); %WINDOWED SIGNAL
            V=abs(WW(k,:)*w_signal.');
            [v(2),indexx(2)]=max(V);
            
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(frac_wind_index+1,:)); %WINDOWED SIGNAL
            V=abs(WW(k,:)*w_signal.');
            [v(3),indexx(3)]=max(V);
            [v_m,ind]=max(v);
            f0=k(indexx(ind));%+f0-delta;
            frac_wind_index=frac_wind_index+ind-2;
            
            
            if v_m<Thr*v_oldd
                
                break;
            end
            %v_old=v_m;
            
            IF(t_index)=f0;
        end
    end
    IF=IF/(length(Sig));
    %figure; plot(IF)
    % [s1,~,~] = ICCD(Sig,1,IF,orderIF,orderamp,alpha);
    
    % Code to extract the compoennt
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    Sigg=Sig.*s_dechirp;
    
    s2=fftshift(fft(Sigg)); 
    s1=zeros(1,length(s2));
   
    s1(length(Sig)/2-1*delta:length(Sig)/2+1*delta)=s2(length(Sig)/2-1*delta:length(Sig)/2+1*delta);%=1;
    s1=ifft(ifftshift(s1));%.*conj(s_dechirp);
    s1a=s2;
    s1a(length(Sig)/2-1*delta:length(Sig)/2+1*delta)=0;
    s1a=ifft(ifftshift(s1a));%.*conj(s_dechirp);
    
%    s11(length(Sig)/2-11*delta:length(Sig)/2+11*delta)=s2(length(Sig)/2-11*delta:length(Sig)/2+11*delta);%=1;
 %       s11=hamming(length(s2)).'.*s11;

    %s11=hamming(length(s2)).'.*s2;
  %  s11=ifft(ifftshift(s11));%.*conj(s_dechirp);
        s11=Sigg;
   
    s11(indexs)=recover_component(s11(iiii).',iiii,indexs).';
   s2=fftshift(fft(s11)); 
    s12=zeros(1,length(s2));
    s12(length(Sig)/2-2*delta:length(Sig)/2+2*delta)=s2(length(Sig)/2-2*delta:length(Sig)/2+2*delta);%.*hamming(4*delta+1).';%=1;
    s11=ifft(ifftshift(s12));%.*conj(s_dechirp);
    
    
    
    %s1=s1.*conj(s_dechirp);
    s11=s11.*conj(s_dechirp);
    
    Sig(iiii)=Sig(iiii)-s11(iiii);
    sig_out=sig_out+s11;
    %Sig=s1a.*conj(s_dechirp);
    %Sig(indexs)=0;
    if sum(abs(s11.^2))<e_max*thr
        break;
    else
        %  e_max
        if sum(abs(s11.^2))>e_max
            e_max=sum(abs(s11.^2));
        end
        fidexmult(iii,:) = IF;
    end
    
end
sig_out(iiii)=Siggg(iiii);

end


