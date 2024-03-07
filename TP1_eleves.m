
% Script for computing the BER for BPSK/QPSK modulation in ISI Channels
% 
close all;
clear all;

%% Simulation parameters
% On décrit ci-après les paramètres généraux de la simulation

%Frame length
M=4; %2:BPSK, 4: QPSK
N  = 10000; % Number of transmitted bits or symbols
Es_N0_dB = [0:30]; % Eb/N0 values
%Multipath channel parameters
hc=[1 0.8*exp(1i*pi/3) 0.3*exp(1i*pi/6) ];%0.1*exp(1i*pi/12)];%ISI channel
%hc=[0.04, -0.05, 0.07, -0.21, -0.5, 0.72, 0.36, 0, 0.21, 0.03, 0.07];
% a=0.8;
% hc=[1 -a];
Lc=length(hc);%Channel length
ChannelDelay=0; %delay is equal to number of non causal taps

%V i t e r b i d e c o d i n g p a r a m e t e r s
const = qammod(( 0 : M - 1 )' , M)/sqrt(2); %r e f e r e n c e Gray QPSK c o n s t e l l a t i o n
tblen = 10 ; %T r a c e b a c k d e p t h
nsamp = 1 ; %O v e r s a m p l i n g r a t e
preamble = [ ] ;
postamble = [ ] ;

%Preallocations
nErr_zfinf=zeros(1,length(Es_N0_dB));
for ii = 1:length(Es_N0_dB)

   % BPSK symbol generations
%    bits = rand(1,N)>0.5; % generating 0,1 with equal probability
%    s = 1-2*bits; % BPSK modulation following: {0 -> +1; 1 -> -1} 
   
    % QPSK symbol generations
   bits = rand(2,N)>0.5; % generating 0,1 with equal probability
   s = 1/sqrt(2)*((1-2*bits(1,:))+1j*(1-2*bits(2,:))); % QPSK modulation following the BPSK rule for each quadatrure component: {0 -> +1; 1 -> -1} 
   sigs2=var(s);
   
   % Channel convolution: equivalent symbol based representation
   z = conv(hc,s);  
   
   %Generating noise
   sig2b=10^(-Es_N0_dB(ii)/10);
   %n = sqrt(sig2b)*randn(1,N+Lc-1); % white gaussian noise, BPSK Case
    n = sqrt(sig2b/2)*randn(1,N+Lc-1)+1j*sqrt(sig2b/2)*randn(1,N+Lc-1); % white gaussian noise, QPSK case
   
   % Adding Noise
   y = z + n; % additive white gaussian noise

   %% zero forcing equalization
   % We now study ZF equalization
   
   %Unconstrained ZF equalization, only if stable inverse filtering
   
   
   %%
   % 
   %  The unconstrained ZF equalizer, when existing is given by 
   % 
   % $w_{,\infty,zf}=\frac{1}{h(z)}$
   % 
   %%
   % 
   
    s_zf=filter(1,hc,y);%if stable causal filter is existing
    bhat_zf = zeros(2,length(bits));
    bhat_zf(1,:)= real(s_zf(1:N)) < 0;
    bhat_zf(2,:)= imag(s_zf(1:N)) < 0;
    nErr_zfinfdirectimp(1,ii) = size(find([bits(:)- bhat_zf(:)]),1);
    %Otherwise, to handle the non causal case
    Nzf=200;
    [r, p, k]=residuez(1, hc);
    [w_zfinf]=ComputeRI( Nzf, r, p, k );
    s_zf=conv(w_zfinf,y);
    bhat_zf = zeros(2,length(bits));
    bhat_zf(1,:)= real(s_zf(Nzf:N+Nzf-1)) < 0;
    bhat_zf(2,:)= imag(s_zf(Nzf:N+Nzf-1)) < 0;
    
    nErr_zfinf(1,ii) = size(find([bits(:)- bhat_zf(:)]),1);


    %% MMSE equalization
    
    deltac= zeros( 1 , 2 * Lc - 1 ) ;
    deltac ( Lc ) = 1 ;
    Nmmse = 200 ;%c a u s a l p a r t
    [ r , p , k ] = residuez ( fliplr ( conj ( hc ) ) , ( conv ( hc , fliplr ( conj ( hc ) ) ) + ( sig2b / sigs2) * deltac ) ) ;
    [ w_mmseinf ] = ComputeRI ( Nmmse , r , p , k ) ;
    s_mmse=conv(w_mmseinf,y);
    bhat_mmse = zeros(2,length(bits));
    bhat_mmse(1,:)= real(s_mmse(Nmmse:N+Nmmse-1)) < 0;
    bhat_mmse(2,:)= imag(s_mmse(Nmmse:N+Nmmse-1)) < 0;
    
    nErr_mmseinf(1,ii) = size(find([bits(:)- bhat_mmse(:)]),1);

    %% ZF Contraint
    Nw = 30; % C'est à optimiser 
    H= toeplitz ( [hc(1) zeros(1 , Nw -1)]' , [hc , zeros(1 , Nw-1)]);
    Ry = (conj(H)*H.');
    p = zeros(Nw+Lc - 1, 1);

    P = (H.'*inv(Ry))*conj(H);
    [alpha, dopt] = max(diag(abs(P)));

    p(dopt) = 1;
    Gamma = conj(H)*p;

    w_zf_ls = (inv(Ry)*Gamma).';

    sig_ei_opt = sigs2 - conj(w_zf_ls)*Gamma;
    bias = 1 - sig_ei_opt/sigs2;
    shat = conv(w_zf_ls, y);
    shat = shat(dopt:end);

    bhat_zf_ls = zeros(2,length(bits));
    bhat_zf_ls(1,:)= real(shat(1:N)) < 0;
    bhat_zf_ls(2,:)= imag(shat(1:N)) < 0;

    nErr_zf_lsinf(1,ii) = size(find([bits(:)- bhat_zf_ls(:)]),1);

    %% MMSE Contraint

    Nw = 100;
    H= toeplitz ( [hc(1) zeros(1 , Nw -1)]' , [hc , zeros(1 , Nw-1)]);
    Ry = sigs2*(conj(H)*H.') + sig2b*eye(Nw);
    p = zeros(Nw+Lc - 1, 1);

    P = (H.'*inv((Ry/sigs2))*conj(H));
    [alpha, dopt] = max(diag(abs(P)));

    p(dopt) = 1;
    Gamma = conj(H)*p;

    w_mmse_ls = (inv(Ry)*Gamma).';

    sig_ei_opt = sigs2 - conj(w_mmse_ls)*Gamma;
    bias = 1 - sig_ei_opt/sigs2;
    shat = conv(w_mmse_ls, y);
    shat = shat(dopt:end);

    bhat_mmse_ls = zeros(2,length(bits));
    bhat_mmse_ls(1,:)= real(shat(1:N)) < 0;
    bhat_mmse_ls(2,:)= imag(shat(1:N)) < 0;

    nErr_mmse_lsinf(1,ii) = size(find([bits(:)- bhat_mmse_ls(:)]),1);

    %% ML
    s_ml = mlseeq( y, hc, const , tblen , 'rst' , nsamp , [] , []);

    bhat_ml = zeros(2,length(bits));
    bhat_ml(1,:)= real(s_ml(1:N)) < 0;
    bhat_ml(2,:)= imag(s_ml(1:N)) < 0;

    nErr_ml(1,ii) = size(find([bits(:)- bhat_ml(:)]),1);

end

simBer_zfinfdirectimp = nErr_zfinfdirectimp/N/log2(M); % simulated ber

simBer_zfinf = nErr_zfinf/N/log2(M); % simulated ber

simBer_mmseinf = nErr_mmseinf/N/log2(M); 

simBer_zf_lsinf = nErr_zf_lsinf/N/log2(M);

simBer_mmse_lsinf = nErr_mmse_lsinf/N/log2(M);

simBer_ml = nErr_ml/N/log2(M);

% plot

figure
semilogy(Es_N0_dB,simBer_zfinfdirectimp(1,:),'bs-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_zfinf(1,:),'rs-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_mmseinf(1,:),'gs-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_zf_lsinf(1,:),'ys-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_mmse_lsinf(1,:),'ks-','Linewidth',2);
hold on
semilogy(Es_N0_dB,simBer_ml(1,:),'cs-','Linewidth',2);
axis([0 50 10^-6 0.5])
grid on
legend('sim-zf-inf/direct','sim-zf-inf/direct', 'sim-mmse-inf', 'sim-zf-ls', 'sim-mmse-ls', 'sim-ml');
xlabel('E_s/N_0, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for QPSK in ISI with ZF equalizers')

figure
title('Impulse response')
stem(real(w_zfinf))
hold on
stem(real(w_zfinf),'r-')
ylabel('Amplitude');
xlabel('time index')


figure
title('Impulse response')
stem(real(w_mmseinf))
hold on
stem(real(w_mmseinf),'r-')
ylabel('Amplitude');
xlabel('time index')

figure
title('Impulse response')
stem(real(w_zf_ls))
hold on
stem(real(w_zf_ls),'r-')
ylabel('Amplitude');
xlabel('time index')

figure
title('Impulse response')
stem(real(w_mmse_ls))
hold on
stem(real(w_mmse_ls),'r-')
ylabel('Amplitude');
xlabel('time index')
