% Auteurs : EL MARZOUGUI OMAR & HAKIM BASSBASSE
% Fichier Matlab : Transmisstion multi-utilisateurs : SC-FDMA

close all;
clear;

%% Paramètres de simulation

% Paramètres de la modulations
M = 4; %Modulation order 
Nframe = 100;
Nfft=1024;
Ncp=8;
Ns=Nframe*(Nfft+Ncp);
N= log2(M)*Nframe*Nfft;

M_subcarrier = 256;

% Paramètres des canaux
Eb_N0_dB = 0:30; % Eb/N0 values

h1=[1 -0.9];
L1=length(h1);%Channel length
H1=fft(h1, Nfft);

h2=[1 0.8*exp(1i*pi/3) 0.3*exp(1i*pi/6) ];
L2 = length(h2);
H2 = fft(h2, Nfft);

n1Err_mmsefde=zeros(1,length(Eb_N0_dB));
n2Err_mmsefde=zeros(1,length(Eb_N0_dB));

%% Chaine de transmission

for ii = 1:length(Eb_N0_dB)

    % Message User 1
    bits1= randi([0 1],N,1);
    s1 = qammod(bits1,M,'InputType','bit');
    sigs12=var(s1);

    % Message User 2
    bits2= randi([0 1],N,1);
    s2 = qammod(bits2,M,'InputType','bit');
    sigs22=var(s2);

    % Serial to Parallel & M-Point DFT
    s1mat = reshape(s1, M_subcarrier, []);
    s1mat = fft(s1mat, M_subcarrier);

    s2mat = reshape(s2, M_subcarrier, []);
    s2mat = fft(s2mat, M_subcarrier);

    % Subcarrier Mapping : Localized
    s1mat_Nfft = zeros(Nfft, length(s1mat(1, 1:end)));
    s1mat_Nfft(30:30+(M_subcarrier) -1, 1:end) = s1mat;

    s2mat_Nfft = zeros(Nfft, length(s1mat(2, 1:end)));
    s2mat_Nfft(300:300+(M_subcarrier) -1, 1:end) = s2mat;

    % N-point ifft
    s1_ifft = ifft(s1mat_Nfft, Nfft);
    s2_ifft = ifft(s2mat_Nfft, Nfft);

    % Prefixe Cyclique
    s1matcp=[s1_ifft(end-Ncp+1:end,:);s1_ifft];
    s1cp=reshape(s1matcp,1,[]); % Signal USER 1

    s2matcp=[s2_ifft(end-Ncp+1:end,:);s2_ifft];
    s2cp=reshape(s2matcp,1,[]); % Signal USER 2

    % Canal user 1 + bruit
    z1 = filter(h1,1,s1cp);  

    sig2b=10^(-Eb_N0_dB(ii)/10);
    n1 = sqrt(sig2b/2)*randn(1,length(z1))+1j*sqrt(sig2b/2)*randn(1,length(z1)); % white gaussian noise, 
   
    y1cp = z1 + n1; % additive white gaussian noise

    % Canal user 2 + bruit 
    z2 = filter(h2,1,s2cp);  

    sig2b=10^(-Eb_N0_dB(ii)/10);
    n2 = sqrt(sig2b/2)*randn(1,length(z2))+1j*sqrt(sig2b/2)*randn(1,length(z2)); % white gaussian noise, 
   
    y2cp = z2 + n2; % additive white gaussian noise   
    
    % Signal à la récéption :
    ycp = y1cp + y2cp;
    
    % Suppression prefixe cyclique
    ycp_mat = reshape(ycp, Ncp +Nfft, []);
    ycp_mat = ycp_mat(Ncp +1:end, :);

    % N-point DFT
    y_fft = fft(ycp_mat, Nfft);
    
    % Egaliseur MMSE  
    W1_mmse = conj(H1)./(abs(H1).^2 + (sig2b/sigs12)); % Egaliseur signal 1
    W2_mmse = conj(H2)./(abs(H2).^2 + (sig2b/sigs22)); % Egaliseur signal 2

    % Egalisation
    S1_mmse = diag(W1_mmse)*y_fft;
    S2_mmse = diag(W2_mmse)*y_fft;

    % Subcarrier Demapping
    S1_mmse = S1_mmse(30:30+(M_subcarrier) -1, 1:end);
    S2_mmse = S2_mmse(300:300+(M_subcarrier) -1, 1:end);

    % M-point IDFT
    s1_mmse = ifft(S1_mmse, M_subcarrier);
    s2_mmse = ifft(S2_mmse, M_subcarrier);

    %Detection
    bits1_mmseeq = qamdemod(s1_mmse, M, 'OutputType','bit');
    n1Err_mmsefde(1,ii) = size(find(bits1(:)- bits1_mmseeq(:)),1);

    bits2_mmseeq = qamdemod(s2_mmse, M, 'OutputType','bit');
    n2Err_mmsefde(1,ii) = size(find(bits2(:)- bits2_mmseeq(:)),1);
end

%% Performances

simBer_s1 = n1Err_mmsefde/N; % simulated ber
simBer_s2 = n2Err_mmsefde/N; % simulated ber

% plot

figure
semilogy(Eb_N0_dB,simBer_s1(1,:),'bs-','Linewidth',2);
hold on
semilogy(Eb_N0_dB,simBer_s2(1,:),'rd-','Linewidth',2);
axis([0 50 10^-6 0.5])
grid on
legend('sim-s1','sim-s2');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Performances en transmission multi-utilisateurs : SC-FDMA')