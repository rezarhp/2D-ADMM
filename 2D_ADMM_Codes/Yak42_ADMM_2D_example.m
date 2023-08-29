%---------------------------------------------------------------------------
% Code Implementation by Hamid Reza Hashempour
% This is a simple sample code related to the following paper:
% "Sparsity-Driven ISAR Imaging Based on Two-Dimensional ADMM"
% H. R. Hashempour, IEEE Sensors Journal, vol. 20, no. 22, pp. 13349-13356, 15 Nov. 2020.
% DOI: 10.1109/JSEN.2020.3006105
% Please kindly cite the paper when utilizing this code in your work.
% This code provides implementations for the 2D-SL0, GP-SOONE, 2D-ADMM, and Fast 2D-ADMM algorithms.
% Note: This code is provided for illustrative purposes. It has not been optimized and the output is not guaranteed.
%---------------------------------------------------------------------------

clc;
clear;
close all;
%%
load Yak42;

prf=100;

y=y(:,64*2+1:64*3);
y    = y./max(max(y));

[N,M] = size(y);
[na,nr] = size(y);     %
nstd    = .505*0.03*sqrt(2);    %1.57 
noise = random('normal',0,nstd,na,nr) + 1j*random('normal',0,nstd,na,nr) ;
SNR = db(norm(y,'fro')^2/(var(noise(:))*na*nr))/2;
fprintf('SNR=%d;\n',SNR)
y=y+noise;
SNR=round(SNR)
sig = ifft(y);
y=fft(sig,512);

n=0:N-1;
P=2*N;
k=(0:P-1).';
Fr=1/sqrt(P)*exp(1i*2*pi*k*n/P).';
m=0:M-1;
Q=2*M;%%%%%%%%%%%%%%%%%%%%%%%%%%2*M;
k=(0:Q-1).';
Fa=1/sqrt(Q)*exp(-1i*2*pi*k*m/Q).';

fd   = linspace(-prf/2,prf/2,Q);
raxis = linspace(-N/2*0.375,N/2*0.375,P);

% tic
% im1=Fr'*sig*conj(Fa);
% toc
figure
contour(1:64,1:256,abs(sig),40);
xlabel('Azimuth samples');ylabel('Range-freq. samples');
set(gca,'FontName', 'Arial', 'FontSize',14);

tic
im2 = pinv(Fr)*sig*pinv(Fa).';
toc
figure
contour(fd,raxis,abs(fftshift(im2,2)),40);
set(gca,'FontName', 'Arial', 'FontSize',14);
xlabel('Doppler (Hz) ');ylabel('Range (m)');
% title(['P=', num2str(Q)])
xlim([-40 40])
ylim([-35 35])


IE(1) = Entropy_img(im2);


tic
s3=SL0_2D(Fr,Fa,sig,.007,.5,2,3);%.005
time_2D_SL0=toc
IE(2) = Entropy_img(s3);


figure
contour(fd,raxis,abs(fftshift(s3,2)),40);
set(gca,'FontName', 'Arial', 'FontSize',14);
xlabel('Doppler (Hz) ');ylabel('Range (m)');
% title(['P=', num2str(Q)])
xlim([-40 40])
ylim([-35 35])
% saveas(gcf,['sl02D', num2str(SNR)],'epsc')
% saveas(gcf,['sl02D', num2str(SNR),'.fig'])

L0=10;L1=8;L=10;sigma_decrease_factor=.5;sigma_min=.048;
tic
s=GP_SOONE(Fr,Fa,sig, sigma_min, sigma_decrease_factor,L0,L1, L);
time_GP_SOONE=toc

figure
contour(fd,raxis,abs(fftshift(s,2)),40);
set(gca,'FontName', 'Arial', 'FontSize',14);
xlabel('Doppler (Hz) ');ylabel('Range (m)');
% title(['P=', num2str(Q)])
xlim([-40 40])
ylim([-35 35])
% saveas(gcf,['soone2D', num2str(SNR)],'epsc')
% saveas(gcf,['soone2D', num2str(SNR),'.fig'])

IE(3) = Entropy_img(s)


%==========================================================================
%%
error=1e-5;
alpha=.00650;%.0045;55;.008;%0.0055;
tic
im=admm_2D_fast(Fr,Fa,sig,error,alpha);
time_admm_2D_fast=toc
IE(4) = Entropy_img( im+eps );

tic
im=admm_2D(Fr,Fa,sig,error,alpha);
time_admm_2D=toc
IE(5) = Entropy_img( im+eps )

figure
contour(fd,raxis,abs(fftshift(im,2)),40);
set(gca,'FontName', 'Arial', 'FontSize',14);
xlabel('Doppler (Hz) ');ylabel('Range (m)');
% title(['P=', num2str(Q)])
xlim([-40 40])
ylim([-35 35])


