% Script for computing the BER for QAM modulation in ISI Channels with FDE
% 
close all;
clear all;

%% Simulation parameters
% On d�crit ci-apr�s les param�tres g�n�raux de la simulation

%modulation parameters
M = 4; %Modulation order 
Nframe = 100;
Nfft=1024;
Ncp=8;
Ns=Nframe*(Nfft+Ncp);
N= log2(M)*Nframe*Nfft;

%Channel Parameters
Eb_N0_dB = [0:30]; % Eb/N0 values
%Multipath channel parameters
hc=[1 -0.9];
Lc=length(hc);%Channel length
H=fft(hc,Nfft);
%Preallocations
nErr_zffde=zeros(1,length(Eb_N0_dB));
nErr_mmsefde=zeros(1,length(Eb_N0_dB));
for ii = 1:length(Eb_N0_dB)

   %Message generation
   bits= randi([0 1],N,1);
   s = qammod(bits,M,'InputType','bit');
   sigs2=var(s);
   
   %Add CP
   smat=reshape(s,Nfft,Nframe);
   smatcp=[smat(end-Ncp+1:end,:);smat];
   scp=reshape(smatcp,1,(Nfft+Ncp)*Nframe);
   
    % Channel convolution: equivalent symbol based representation
   z = filter(hc,1,scp);  
   
   %Generating noise
   sig2b=10^(-Eb_N0_dB(ii)/10);
   
   %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, 
   n = sqrt(sig2b/2)*randn(1,Ns)+1j*sqrt(sig2b/2)*randn(1,Ns); % white gaussian noise, 
   
    % Noise addition
   ycp = z + n; % additive white gaussian noise

   %remove CP
   ycp_mat = reshape(ycp, Ncp +Nfft, Nframe);
   ycp_mat = ycp_mat(Ncp +1:end, :);
   
   %% ZF 
   %FDE
   W_zf = 1./H;

   %DFT
   ycp_fft = fft(ycp_mat,Nfft);

   %Egalisation
   S_zf = diag(W_zf)*ycp_fft;
   s_zf = ifft(S_zf, Nfft);

   %Detection
   bhat_zfeq = qamdemod(s_zf, M, 'OutputType','bit');
   nErr_zffde(1,ii) = size(find([bits(:)- bhat_zfeq(:)]),1);
   

   %% MMSE  
   W_mmse = conj(H)./(abs(H).^2 + (sig2b/sigs2));

   %Egalisation
   S_mmse = diag(W_mmse)*ycp_fft;
   s_mmse = ifft(S_mmse, Nfft);

   %Detection
   bhat_mmseeq = qamdemod(s_mmse, M, 'OutputType','bit');
   nErr_mmsefde(1,ii) = size(find([bits(:)- bhat_mmseeq(:)]),1);

end

simBer_zf = nErr_zffde/N; % simulated ber
simBer_mmse = nErr_mmsefde/N; % simulated ber

% plot

figure
semilogy(Eb_N0_dB,simBer_zf(1,:),'bs-','Linewidth',2);
hold on
semilogy(Eb_N0_dB,simBer_mmse(1,:),'rd-','Linewidth',2);
axis([0 70 10^-6 0.5])
grid on
legend('sim-zf-fde','sim-mmse-fde');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK in ISI with ZF and MMSE equalizers')


