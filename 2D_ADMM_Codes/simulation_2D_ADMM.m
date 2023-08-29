%----------------------------------------------------------------
% Analysis of SGP-SOONE vs SNR for a synthetic target
%----------------------------------------------------------------
%---------------------------------------------------------------------------
% Code Implementation by Hamid Reza Hashempour
% This is a simple sample code related to the following paper:
% "Sparsity-Driven ISAR Imaging Based on Two-Dimensional ADMM"
% H. R. Hashempour, IEEE Sensors Journal, vol. 20, no. 22, pp. 13349-13356, 15 Nov. 2020.
% DOI: 10.1109/JSEN.2020.3006105
% Please kindly cite the paper when utilizing this code in your work.
% This code provides implementations for the 2D-SL0, GP-SOONE, 2D-ADMM, and algorithms.
% Note: This code is provided for illustrative purposes. It has not been optimized and the output is not guaranteed.
%---------------------------------------------------------------------------



clc
clear
close all
%%
iter=1;200;
%---Radar parameters------------------------------------------------
c = 3e8;               % speed of EM wave [m/s]
fc=10e9;               % Center frequency
lambda=c/fc;
T1 = 2e-6;              % duration of single batch [s]
Tp = T1;
PRF=50;
T2 =1/PRF;              % duration of PRI [s]
BW=500e6;
fs=500e6;
N=50;
M=50;
SNR_db = 10;-10:5:30;
noise=1;

MSE_FFT_2D1=zeros(iter,length(SNR_db));
MSE_GP_SOONE1=zeros(iter,length(SNR_db));
MSE_SL0_2D1=zeros(iter,length(SNR_db));
MSE_ADMM_2D1=zeros(iter,length(SNR_db));

for idx=1:iter
    for index1=1:length(SNR_db)
        %---target parameters------------------------------------------------
        W = .05;%%%%%% 0.0390625 for N=64            % Angular velocity [rad/s]
        Vr1 =0;100;500;10;        % radial translational  velocity of EM wave [m/s]
        ar1 =0;10;.05;            % radial accelation of EM wave [m/s^2]
        R0 = 5e3;             % target's initial distance from radar [m]
        theta0 = 0;           % Look angle of the target
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cr_res= lambda/(2*abs(W)*N*T2);
        r_res=c/(2*fs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K= M;
        L= N;
        target=zeros(K,L);
        target(round(K/2),round([L/2-16 L/2-8 L/2 L/2+8 L/2+16]))=255;
        target(round([K/2-6 K/2+6]),round(L/2-12))=255;
        target(round([K/2-16 K/2-8 K/2+8 K/2+16]),round(L/2+12))=255;
        
        
        
        [K, L]=size(target); ntarget=K*L;
        tnum=1; xn=zeros(ntarget,1); yn=xn; Fn=xn; % Target Intialization Variables
        for m=1:K
            for n=1:L
                xn(tnum)=-(n-L/2);
                yn(tnum)=(m-K/2);
                Fn(tnum)=double(target(m,n))/255;
                tnum=tnum+1;
            end
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Xc=Xc;
        % Yc=Yc;
        %         target_scene=reshape(Fn,K,L);
%         craxis   = linspace(-N/2*cr_res,N/2*cr_res,N);
%         raxis    = linspace(-M/2*r_res,M/2*r_res,M);
%         figure;imagesc(craxis,raxis,abs(reshape(Fn,K,L)))
%         title('Original scene')
%         set(gca,'FontName', 'Arial', 'FontSize',14);
%         ylabel('Cross-range (m) ');xlabel('Range (m)');
%         saveas(gcf,'orig_scene','epsc')
%         saveas(gcf,['orig_scene','.fig'])
        %---Figure --------------------------------------------------------
        Xc=xn(Fn~=0)*r_res;
        Yc=yn(Fn~=0)*r_res;
%         h = figure;
%         plot(Yc,Xc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0]);
%         % set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
%         ggg=max(max(abs(Yc),abs(Xc)))+5;
%         axis([-ggg ggg -ggg ggg])
%         xlabel('X [m]'); ylabel('Y [m]');
        
        T = N*T2;               % Dwell time [s]
        
        
        
        
        %===========================================================
        %Scattering centers in cylindirical coordinates
        [theta,r] = cart2pol(Xc,Yc);
        theta = theta+theta0*0.017455329; %add initial angle
        
        
        n=1:N;
        T11 = 2*R0(1)/c+(n-1)*T2;%T1/2+
        
        Rvr(1,:) = Vr1*T11+(0.5*ar1)*(T11.^2);%+.005*sin(2*pi*2*T11);%Range Displacement due to radial vel. & acc.
        Tetw = W*T11;% Rotational Displacement due to angular vel.
        
        
        %--- sampling & time parameters -----------------------------
        dt = 1/fs;                  % sampling time interval
        Kchirp = BW/T1;              % chirp pulse parameter
        Ts=(2*(R0-(M/2-1)*r_res))/c; % Start time of sampling
        Tf=(2*(R0+M/2*r_res))/c+Tp; % End time of sampling
        %..................Measurement Parameters
        rbins=round(((Tf-Ts))/dt); % Number of time (Range) samples
        t=Ts+(0:rbins-1)*dt; % Time array for data acquisition
        %---------------------------------------------
        tt=(0:dt:Tp-dt)-Tp/2;
        s=exp(1j*pi*Kchirp*tt.^2);
        %---------------------------------------------
        M1 =rbins;%round(T2*fs1);      % range samples
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Es =zeros(N,M1);
        
        TargetPSample=zeros(size(Xc));
        
        tic
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        u=zeros(N,M1);
        % T11=zeros(N,M);
        % for n=1:N
        %     T11(n,:) = 2*R0(1)/c+(n-1)*T2+(0:M-1)*1/fs;%T1/2+
        % end
        % Rvr = Vr1*T11+(0.5*ar1)*(T11.^2);%+.005*sin(2*pi*2*T11);%Range Displacement due to radial vel. & acc.
        % Tetw = W*T11;% Rotational Displacement due to angular vel.
        % kk=0:M1-1;
        for m =1: length(Xc)
            Es =zeros(N,M1);
            %         R = R0 + Rvr + Yc(m)*cos(Tetw) + Xc(m)*sin(Tetw);%
            R = R0 + Rvr + Yc(m) + Xc(m)*Tetw;%*cos(Tetw) *sin(Tetw)
            for n=1:N
                td=t-2*(R0 + Rvr(n) + Yc(m))/c;
                %         td=t-2*(R(n))/c;
                
                Es(n,1:M1) = Es(n,1:M1)+exp(-1j*(4*pi*fc*R(n)/c)+1j*pi*Kchirp*...
                    ((td-Tp/2).^2)).*(td >= 0 & td <= Tp);
            end
            
            u=u+Es;
        end
%         figure;imagesc(abs(u))
        
        if noise==1
            E_signal = sum(sum(abs(u.^2)))/(M1*N);
            SNR=10^(SNR_db(index1)/10);
            std_dev=sqrt(E_signal/SNR);
            s_n=std_dev*randn(size(u));
            u=u+s_n;
        end
        
%         figure;imagesc(abs(u))
        
        
        % u1=u;
        M2=round(Tp/dt);
        u1=zeros(N,M1-M2+1);
        for ii=1:N
            u1(ii,:)=conv(u(ii,:),conj(fliplr(s)),'valid');%.*taylorwin(M,5,-35).'
        end
%         figure;imagesc(db(u1))
%         xlabel('Range, sample')
%         ylabel('azimuth, sample')
%         
%         figure;imagesc(abs(fftshift(fft(u1.*(taylorwin(N,5,-35)*ones(1,M1-M2+1))),1)))
%         xlabel('Range, sample')
%         ylabel('azimuth, sample')
        
        
        %==========================================================================
        P=2*N;
        n=0:N-1;
        k=(0:P-1).';
        F=1/sqrt(P)*exp(-1i*2*pi*k*n/P).';
        
        Q=2*M;
        m=0:M-1;
        k=(0:Q-1).';
        F1=1/sqrt(Q)*exp(-1i*2*pi*k*m/Q).';
        L0=5;L1=8;L=3;sigma_decrease_factor=.5;sigma_min=.1;
        x1=fft(u1,[],2);
        mu=2;
        % s0 = pinv(F)*x1*pinv(F1).';
        s0 = (F')*x1*conj(F1);
        s1=GP_SOONE(F,F1,x1, sigma_min, sigma_decrease_factor,L0,L1, L);
        s2=SL0_2D(F,F1,x1,sigma_min, sigma_decrease_factor,mu,L);%*1000000/10
        error=1e-4;
        alpha=.5;.008;%0.0055;
        x1=x1./max(abs(x1(:)));
        s3=admm_2D(F,F1,x1,error,alpha);
        % A=kron(F1,F);
        % y=x1(:);
        % s3=SL0(A,y,sigma_min*1000000/10, sigma_decrease_factor,mu,L);
        %%
        fd   = linspace(-PRF/2,PRF/2,P);
        raxis = linspace(-Q/2*r_res/2,Q/2*r_res/2,Q);
        s0_norm=abs(fftshift(s0,1))/max(abs(s0(:)));
        figure;imagesc(raxis,fd,(abs(s0_norm)))
        title('2D-FFT')
        set(gca,'FontName', 'Arial', 'FontSize',14);
        ylabel('Doppler (Hz) ');xlabel('Range (m)');
        axis xy
%         saveas(gcf,['fftsim', num2str(SNR_db(idx))],'epsc')
%         saveas(gcf,['fftsim', num2str(SNR_db(idx)),'.fig'])

        s1_norm=abs(fftshift(s1,1))/max(abs(s1(:)));
        figure;imagesc(raxis,fd,s1_norm)
        title('2D-GP-SOONE')
        set(gca,'FontName', 'Arial', 'FontSize',14);
        ylabel('Doppler (Hz) ');xlabel('Range (m)');
        axis xy
%         saveas(gcf,['soonesim', num2str(SNR_db(idx))],'epsc')
%         saveas(gcf,['soonesim', num2str(SNR_db(idx)),'.fig'])
        
        s2_norm=abs(fftshift(s2,1))/max(abs(s2(:)));
        figure;imagesc(raxis,fd,s2_norm)
        title('2D-SL0')
        set(gca,'FontName', 'Arial', 'FontSize',14);
        ylabel('Doppler (Hz) ');xlabel('Range (m)');
        axis xy
%         saveas(gcf,['sl0sim', num2str(SNR_db(idx))],'epsc')
%         saveas(gcf,['sl0sim', num2str(SNR_db(idx)),'.fig'])
                
        s3_norm=abs(fftshift(s3,1))/max(abs(s3(:)));
        figure;imagesc(raxis,fd,s3_norm)
        title('2D-ADMM')
        set(gca,'FontName', 'Arial', 'FontSize',14);
        ylabel('Doppler (Hz) ');xlabel('Range (m)');
        axis xy
%         saveas(gcf,['admmsim', num2str(SNR_db(idx))],'epsc')
%         saveas(gcf,['admmsim', num2str(SNR_db(idx)),'.fig'])
        
        % figure;imagesc(abs(fftshift(reshape(s3,100,100),1)))
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K= P;
        L= Q;
        target=zeros(K,L);
        center=[K/2-1,L/2+1];
        target(center(1),round([center(2)-16*2 center(2)-8*2 center(2) center(2)+8*2 center(2)+16*2]))=255;
        target(round([center(1)-6*2 center(1)+6*2]),round(center(2)+12*2))=255;
        target(round([center(1)-16*2 center(1)-8*2 center(1)+8*2 center(1)+16*2]),round(center(2)-12*2))=255;
        
        [K, L]=size(target); ntarget=K*L;
        tnum=1; xn=zeros(ntarget,1); yn=xn; Fn=xn; % Target Intialization Variables
        for m=1:K
            for n=1:L
                xn(tnum)=-(n-L/2);
                yn(tnum)=(m-K/2);
                Fn(tnum)=double(target(m,n))/255;
                tnum=tnum+1;
            end
        end
        target_scene=reshape(Fn,K,L);
%         figure;imagesc(abs(target_scene))
%         title('Original Target scene')
         craxis   = linspace(-P/2*cr_res/2,P/2*cr_res/2,P);
        raxis    = linspace(-Q/2*r_res/2,Q/2*r_res/2,Q);
        figure;imagesc(craxis,raxis,abs(target_scene))
        title('Original scene')
        set(gca,'FontName', 'Arial', 'FontSize',14);
        ylabel('Cross-range (m) ');xlabel('Range (m)');
        axis xy
%         saveas(gcf,'orig_scene','epsc')
%         saveas(gcf,['orig_scene','.fig'])
        % MSE_FFT_2D(index1)=(norm(s0_norm(:)-target_scene(:)));
        % MSE_GP_SOONE(index1)=(norm(s1_norm(:)-target_scene(:)));
        % MSE_SL0_2D(index1)=(norm(s2_norm(:)-target_scene(:)));
        % MSE_ADMM_2D(index1)=(norm(s3_norm(:)-target_scene(:)));
        
        MSE_FFT_2D1(idx,index1)=10*log10(norm(s0_norm(:)-target_scene(:)).^2);
        MSE_GP_SOONE1(idx,index1)=10*log10((norm(s1_norm(:)-target_scene(:))).^2);
        MSE_SL0_2D1(idx,index1)=10*log10((norm(s2_norm(:)-target_scene(:))).^2);
        MSE_ADMM_2D1(idx,index1)=10*log10((norm(s3_norm(:)-target_scene(:))).^2);
    end
    idx
end
%%
% figure;semilogy(undersamp_ratio,mean(MSE_GP_SOONE).^2/(P*Q))
% hold on;semilogy(undersamp_ratio,mean(MSE_SL0_2D).^2/(P*Q))
% hold on;semilogy(undersamp_ratio,mean(MSE_ADMM_2D).^2/(P*Q),'o')
% hold on;semilogy(undersamp_ratio,mean(MSE_FFT_2D).^2/(P*Q))

figure;plot(SNR_db,mean(MSE_FFT_2D1),'-pentagram')
hold on;plot(SNR_db,mean(MSE_GP_SOONE1),'-*')
hold on;plot(SNR_db,mean(MSE_SL0_2D1),'-o')
hold on;plot(SNR_db,mean(MSE_ADMM_2D1),'-square')
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold');
grid
xlabel('SNR [dB]')
ylabel('NMSE [dB]')
legend('2D-FFT','2D-GP-SOONE','2D-SL0','2D-ADMM')

