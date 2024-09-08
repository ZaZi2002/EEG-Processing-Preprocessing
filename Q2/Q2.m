clc
close all
clear all

%% Part 1
load('Ex2.mat');
load("Electrodes.mat");
fs = 200;
t = 0.0001:1/fs:51.2;

%%%%% Calculating energies
p_sig = 0;
p_noise1 = 0;
p_noise2 = 0;
for i = 1:length(X_org(:,1))
    for j = 1:length(X_org(1,:))
        p_sig = p_sig + X_org(i,j)^2; % energy of main signal
        p_noise1 = p_noise1 + X_noise_1(i,j)^2; % energy of noise
        p_noise2 = p_noise2 + X_noise_3(i,j)^2; % energy of noise
    end
end

%%%%% Calculating Sigmas
snr = -10;
sigma_snr10_1 = ((p_sig/p_noise1)*10^(snr/-10))^0.5;
sigma_snr10_2 = ((p_sig/p_noise2)*10^(snr/-10))^0.5;
snr = -20;
sigma_snr20_1 = ((p_sig/p_noise1)*10^(snr/-10))^0.5;
sigma_snr20_2 = ((p_sig/p_noise2)*10^(snr/-10))^0.5;

%%%%% Producing noisy signals with given SNRs
X_10_1 = X_org + sigma_snr10_1*X_noise_1; % noisy signal with -10 SNR
X_20_1 = X_org + sigma_snr20_1*X_noise_1; % noisy signal with -20 SNR
X_10_2 = X_org + sigma_snr10_2*X_noise_3; % noisy signal with -10 SNR
X_20_2 = X_org + sigma_snr20_2*X_noise_3; % noisy signal with -20 SNR

%%%%% Plotting noisy signals
%Noise1 with -10 SNR
offset = max(max(abs(X_10_1)))/3 ;
disp_eeg(X_10_1,offset,fs,Electrodes.labels);
title("Signal + Noise 1 with -10 SNR")
xlim('tight');
grid minor

%Noise1 with -20 SNR
offset = max(max(abs(X_20_1)))/3 ;
disp_eeg(X_20_1,offset,fs,Electrodes.labels);
title("Signal + Noise 1 with -20 SNR")
xlim('tight');
grid minor

%Noise3 with -10 SNR
offset = max(max(abs(X_10_2)))/3 ;
disp_eeg(X_10_2,offset,fs,Electrodes.labels);
title("Signal + Noise 3 with -10 SNR")
xlim('tight');
grid minor

%Noise3 with -20 SNR
offset = max(max(abs(X_20_2)))/3 ;
disp_eeg(X_20_2,offset,fs,Electrodes.labels);
title("Signal + Noise 3 with -20 SNR")
xlim('tight');
grid minor

%% Part 2
%%%%% PCA with pca function
ZM_X_10_1 = X_10_1 - mean(X_10_1,2); % Zero mean
ZM_X_20_1 = X_20_1 - mean(X_20_1,2); % Zero mean
ZM_X_10_2 = X_10_2 - mean(X_10_2,2); % Zero mean
ZM_X_20_2 = X_20_2 - mean(X_20_2,2); % Zero mean
[C_10_1,~,L_10_1] = pca(ZM_X_10_1.');
[C_20_1,~,L_20_1] = pca(ZM_X_20_1.');
[C_10_2,~,L_10_2] = pca(ZM_X_10_2.');
[C_20_2,~,L_20_2] = pca(ZM_X_20_2.');

%%%%% Producing sources matrix
D_10_1 = diag((L_10_1.').^(-0.5))* C_10_1.';
D_20_1 = diag((L_20_1.').^(-0.5))* C_20_1.';
D_10_2 = diag((L_10_2.').^(-0.5))* C_10_2.';
D_20_2 = diag((L_20_2.').^(-0.5))* C_20_2.';
pca_sources_10_1 = D_10_1 * ZM_X_10_1;
pca_sources_20_1 = D_20_1 * ZM_X_20_1;
pca_sources_10_2 = D_10_2 * ZM_X_10_2;
pca_sources_20_2 = D_20_2 * ZM_X_20_2;

%%%%% Plotting Sources
%Noise1 with -10 SNR PCA
offset = max(max(abs(pca_sources_10_1)))/3 ;
disp_eeg(pca_sources_10_1,offset,fs);
title("PCA Sources with noise1 (-10 SNR)")
xlim('tight');
grid minor

%Noise1 with -20 SNR PCA
offset = max(max(abs(pca_sources_20_1)))/3 ;
disp_eeg(pca_sources_20_1,offset,fs);
title("PCA Sources with noise1 (-20 SNR)")
xlim('tight');
grid minor

%Noise3 with -10 SNR PCA
offset = max(max(abs(pca_sources_10_2)))/3 ;
disp_eeg(pca_sources_10_2,offset,fs);
title("PCA Sources with noise3 (-10 SNR)")
xlim('tight');
grid minor

%Noise3 with -20 SNR PCA
offset = max(max(abs(pca_sources_20_2)))/3 ;
disp_eeg(pca_sources_20_2,offset,fs);
title("PCA Sources with noise3 (-20 SNR)")
xlim('tight');
grid minor



%%%%% ICA with COM2 
[F_10_1,W_10_1,K_10_1] = COM2R(X_10_1,32);
[F_20_1,W_20_1,K_20_1] = COM2R(X_20_1,32);
[F_10_2,W_10_2,K_10_2] = COM2R(X_10_2,32);
[F_20_2,W_20_2,K_20_2] = COM2R(X_20_2,32);

%%%%% Finding Sources
Sources_10_1 = W_10_1*X_10_1;
Sources_20_1 = W_20_1*X_20_1;
Sources_10_2 = W_10_2*X_10_2;
Sources_20_2 = W_20_2*X_20_2;

%%%%% Plotting Sources
%Noise1 with -10 SNR ICA
offset = max(max(abs(Sources_10_1)))/3 ;
disp_eeg(Sources_10_1,offset,fs);
title("ICA Sources with noise1 (-10 SNR)")
xlim('tight');
grid minor

%Noise1 with -20 SNR ICA
offset = max(max(abs(Sources_20_1)))/3 ;
disp_eeg(Sources_20_1,offset,fs);
title("ICA Sources with noise1 (-20 SNR)")
xlim('tight');
grid minor

%Noise3 with -10 SNR ICA
offset = max(max(abs(Sources_10_2)))/3 ;
disp_eeg(Sources_10_2,offset,fs);
title("ICA Sources with noise3 (-10 SNR)")
xlim('tight');
grid minor

%Noise3 with -20 SNR ICA
offset = max(max(abs(Sources_20_2)))/3 ;
disp_eeg(Sources_20_2,offset,fs);
title("ICA Sources with noise3 (-20 SNR)")
xlim('tight');
grid minor

%% Part 3&4
%%%%% Keeping spiky sources
SelSources_PCA_1 = [2 3 5 16 32];
SelSources_PCA_2 = [1 4 6 10 11 27 32];
SelSources_ICA_1 = [3 6 14 25];
SelSources_ICA_2 = [4 10 12 18 19];

%%%%% Producing denoised signals
inv_D_10_1 = inv(D_10_1);
inv_D_20_1 = inv(D_20_1);
inv_D_10_2 = inv(D_10_2);
inv_D_20_2 = inv(D_20_2);
PCA_X_den_10_1 = inv_D_10_1(:,SelSources_PCA_1)*Sources_10_1(SelSources_PCA_1,:);
PCA_X_den_20_1 = inv_D_20_1(:,SelSources_PCA_1)*Sources_20_1(SelSources_PCA_1,:);
PCA_X_den_10_2 = inv_D_10_2(:,SelSources_PCA_2)*Sources_10_2(SelSources_PCA_2,:);
PCA_X_den_20_2 = inv_D_20_2(:,SelSources_PCA_2)*Sources_20_2(SelSources_PCA_2,:);
ICA_X_den_10_1 = F_10_1(:,SelSources_ICA_1)*Sources_10_1(SelSources_ICA_1,:);
ICA_X_den_20_1 = F_20_1(:,SelSources_ICA_1)*Sources_20_1(SelSources_ICA_1,:);
ICA_X_den_10_2 = F_10_2(:,SelSources_ICA_2)*Sources_10_2(SelSources_ICA_2,:);
ICA_X_den_20_2 = F_20_2(:,SelSources_ICA_2)*Sources_20_2(SelSources_ICA_2,:);

%%%%% Plotting denoised signals
%PCA Denoised1 with -10 SNR
offset = max(max(abs(PCA_X_den_10_1)))/3 ;
disp_eeg(PCA_X_den_10_1,offset,fs,Electrodes.labels);
title("PCA X-den-1 (-10 SNR)")
xlim('tight');
grid minor

%PCA Denoised1 with -20 SNR
offset = max(max(abs(PCA_X_den_20_1)))/3 ;
disp_eeg(PCA_X_den_20_1,offset,fs,Electrodes.labels);
title("PCA X-den-1 (-20 SNR)")
xlim('tight');
grid minor

%PCA Denoised3 with -10 SNR
offset = max(max(abs(PCA_X_den_10_2)))/3 ;
disp_eeg(PCA_X_den_10_2,offset,fs,Electrodes.labels);
title("PCA X-den-2 (-10 SNR)")
xlim('tight');
grid minor

%PCA Denoised3 with -20 SNR
offset = max(max(abs(PCA_X_den_20_2)))/3 ;
disp_eeg(PCA_X_den_20_2,offset,fs,Electrodes.labels);
title("PCA X-den-2 (-20 SNR)")
xlim('tight');
grid minor


%ICA Denoised1 with -10 SNR
offset = max(max(abs(ICA_X_den_10_1)))/3 ;
disp_eeg(ICA_X_den_10_1,offset,fs,Electrodes.labels);
title("ICA X-den-1 (-10 SNR)")
xlim('tight');
grid minor

%ICA Denoised1 with -20 SNR
offset = max(max(abs(ICA_X_den_20_1)))/3 ;
disp_eeg(ICA_X_den_20_1,offset,fs,Electrodes.labels);
title("ICA X-den-1 (-20 SNR)")
xlim('tight');
grid minor

%ICA Denoised3 with -10 SNR
offset = max(max(abs(ICA_X_den_10_2)))/3 ;
disp_eeg(ICA_X_den_10_2,offset,fs,Electrodes.labels);
title("ICA X-den-2 (-10 SNR)")
xlim('tight');
grid minor

%ICA Denoised3 with -20 SNR
offset = max(max(abs(ICA_X_den_20_2)))/3 ;
disp_eeg(ICA_X_den_20_2,offset,fs,Electrodes.labels);
title("ICA X-den-2 (-20 SNR)")
xlim('tight');
grid minor

%% Part 5
% denoised1 -10 SNR
figure("Name","Part5-10-1");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_10_1(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-1 (-10)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,PCA_X_den_10_1(13,:));
grid minor;
xlim('tight');
title('13th channel PCA x-den-1')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,ICA_X_den_10_1(13,:));
grid minor;
xlim('tight');
title('13th channel ICA x-den-1')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_10_1(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-1 (-10)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,PCA_X_den_10_1(24,:));
grid minor;
xlim('tight');
title('24th channel PCA x-den-1')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,ICA_X_den_10_1(24,:));
grid minor;
xlim('tight');
title('24th channel ICA x-den-1')
xlabel('Time (s)')

% denoised1 -20 SNR
figure("Name","Part5-20-1");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_20_1(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-1 (-20)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,PCA_X_den_20_1(13,:));
grid minor;
xlim('tight');
title('13th channel PCA x-den-1')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,ICA_X_den_20_1(13,:));
grid minor;
xlim('tight');
title('13th channel ICA x-den-1')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_20_1(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-1 (-20)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,PCA_X_den_20_1(24,:));
grid minor;
xlim('tight');
title('24th channel PCA x-den-1')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,ICA_X_den_20_1(24,:));
grid minor;
xlim('tight');
title('24th channel ICA x-den-1')
xlabel('Time (s)')

% denoised2 -10 SNR
figure("Name","Part5-10-2");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_10_2(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-2 (-10)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,PCA_X_den_10_2(13,:));
grid minor;
xlim('tight');
title('13th channel PCA x-den-2')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,ICA_X_den_10_2(13,:));
grid minor;
xlim('tight');
title('13th channel ICA x-den-2')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_10_2(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-2 (-10)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,PCA_X_den_10_2(24,:));
grid minor;
xlim('tight');
title('24th channel PCA x-den-2')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,ICA_X_den_10_2(24,:));
grid minor;
xlim('tight');
title('24th channel ICA x-den-2')
xlabel('Time (s)')

% denoised2 -20 SNR
figure("Name","Part5-20-2");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_20_2(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-2 (-20)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,PCA_X_den_20_2(13,:));
grid minor;
xlim('tight');
title('13th channel PCA x-den-2')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,ICA_X_den_20_2(13,:));
grid minor;
xlim('tight');
title('13th channel ICA x-den-2')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_20_2(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-2 (-20)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,PCA_X_den_20_2(24,:));
grid minor;
xlim('tight');
title('24th channel PCA x-den-2')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,ICA_X_den_20_2(24,:));
grid minor;
xlim('tight');
title('24th channel ICA x-den-2')
xlabel('Time (s)')

%% Part 6
PCA_a_10_1 = 0;
PCA_a_10_2 = 0;
PCA_a_20_1 = 0;
PCA_a_20_2 = 0;
ICA_a_10_1 = 0;
ICA_a_10_2 = 0;
ICA_a_20_1 = 0;
ICA_a_20_2 = 0;
b = 0;
for i = 1:length(X_org(:,1))
    for j = 1:length(X_org(1,:))
        PCA_a_10_1 = PCA_a_10_1 + (X_org(i,j)-PCA_X_den_10_1(i,j))^2;
        PCA_a_20_1 = PCA_a_20_1 + (X_org(i,j)-PCA_X_den_20_1(i,j))^2;
        PCA_a_10_2 = PCA_a_10_2 + (X_org(i,j)-PCA_X_den_10_2(i,j))^2;
        PCA_a_20_2 = PCA_a_20_2 + (X_org(i,j)-PCA_X_den_20_2(i,j))^2;

        ICA_a_10_1 = ICA_a_10_1 + (X_org(i,j)-ICA_X_den_10_1(i,j))^2;
        ICA_a_20_1 = ICA_a_20_1 + (X_org(i,j)-ICA_X_den_20_1(i,j))^2;
        ICA_a_10_2 = ICA_a_10_2 + (X_org(i,j)-ICA_X_den_10_2(i,j))^2;
        ICA_a_20_2 = ICA_a_20_2 + (X_org(i,j)-ICA_X_den_20_2(i,j))^2;

        b = b + (X_org(i,j))^2;
    end
end

PCA_RRMSE_10_1 = (PCA_a_10_1/b)^0.5
PCA_RRMSE_20_1 = (PCA_a_20_1/b)^0.5
PCA_RRMSE_10_2 = (PCA_a_10_2/b)^0.5
PCA_RRMSE_20_2 = (PCA_a_20_2/b)^0.5

ICA_RRMSE_10_1 = (ICA_a_10_1/b)^0.5
ICA_RRMSE_20_1 = (ICA_a_20_1/b)^0.5
ICA_RRMSE_10_2 = (ICA_a_10_2/b)^0.5
ICA_RRMSE_20_2 = (ICA_a_20_2/b)^0.5
