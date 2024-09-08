clc
close all
clear all

%% Part 1
signal1 = load('NewData1.mat');
signal2 = load('NewData3.mat');
load("Electrodes.mat")
fs = 250;
t = 0:1/fs:5;

%%%%% Plotting signals
% NewData1
offset = max(max(abs(signal1.EEG_Sig))) ;
disp_eeg(signal1.EEG_Sig,offset,fs,Electrodes.labels);
title("NewData1 EEG Signal")
xlim('tight');
grid minor

% NewData3
offset = max(max(abs(signal2.EEG_Sig))) ;
disp_eeg(signal2.EEG_Sig,offset,fs,Electrodes.labels);
title("NewData3 EEG Signal")
xlim('tight');
grid minor

%% Part 3
%%%%% ICA with COM2 
[F_1,W_1,K_1] = COM2R(signal1.EEG_Sig,32);
[F_2,W_2,K_2] = COM2R(signal2.EEG_Sig,32);

%%%%% Finding Sources
Sources_1 = W_1*signal1.EEG_Sig;
Sources_2 = W_2*signal2.EEG_Sig;

%% Part 4
%%%%% Plotting sources in time domain
% NewData1 Sources
offset = max(max(abs(Sources_1)))/2 ;
disp_eeg(Sources_1,offset,fs);
title("NewData1 Sources")
xlim('tight');
grid minor

% NewData3 Sources
offset = max(max(abs(Sources_2)))/2 ;
disp_eeg(Sources_2,offset,fs);
title("NewData3 Sources")
xlim('tight');
grid minor


%%%%% Plotting sources in frequency domain
figure('Name','Pwelch for each Source of Data1')
for i = 1:7
    for j = 1:3
        n= j + (i-1)*3;
        subplot(7,3,n);
        pwelch(Sources_1(n,:),fs);
        title("Pwelch for Source " + n);
    end
end

figure('Name','Pwelch for each Source of Data2')
for i = 1:7
    for j = 1:3
        n= j + (i-1)*3;
        subplot(7,3,n);
        pwelch(Sources_2(n,:),fs);
        title("Pwelch for Source " + n);
    end
end

figure('Name', 'Pwelch for Sources of Data1');
pwelch(Sources_1.',fs);

figure('Name', 'Pwelch for Sources of Data2');
pwelch(Sources_2.',fs);


%%%%% Plotting topomaps
figure('Name','Topomap for each Source of Data1')
for i = 1:4
    for j = 1:3
        n= j + (i-1)*3;
        subplot(4,3,n);
        plottopomap(Electrodes.X,Electrodes.Y,Electrodes.labels,F_1(:,n))
        title("Topomap for Data1 Source " + n);
    end
end
figure('Name','Topomap for each Source of Data1')
for i = 5:7
    for j = 1:3
        n= j + (i-1)*3 - 12;
        subplot(3,3,n);
        plottopomap(Electrodes.X,Electrodes.Y,Electrodes.labels,F_1(:,n))
        title("Topomap for Data1 Source " + (n+12));
    end
end


figure('Name','Topomap for each Source of Data3')
for i = 1:4
    for j = 1:3
        n= j + (i-1)*3;
        subplot(4,3,n);
        plottopomap(Electrodes.X,Electrodes.Y,Electrodes.labels,F_2(:,n))
        title("Topomap for Data3 Source " + n);
    end
end
figure('Name','Topomap for each Source of Data3')
for i = 5:7
    for j = 1:3
        n= j + (i-1)*3 - 12;
        subplot(3,3,n);
        plottopomap(Electrodes.X,Electrodes.Y,Electrodes.labels,F_2(:,n))
        title("Topomap for Data3 Source " + (n+12));
    end
end

%% Part 5
%%%%% good sources
SelSources_1 = 1:21;
SelSources_1([4, 10]) = []; %Removing artifact sources
SelSources_2 = 1:21;
SelSources_2([1,7,19]) = []; %Removing artifact sources

%%%%% Producing denoised signals
den_signal_1 = F_1(:,SelSources_1)*Sources_1(SelSources_1,:);
den_signal_2 = F_2(:,SelSources_2)*Sources_2(SelSources_2,:);

%%%%% Plotting denoised signals
%ICA Denoised NewData1
offset = max(max(abs(den_signal_1))) ;
disp_eeg(den_signal_1,offset,fs,Electrodes.labels);
title("ICA Denoised NewData1");
xlim('tight');
grid minor

%ICA Denoised NewData2
offset = max(max(abs(den_signal_2))) ;
disp_eeg(den_signal_2,offset,fs,Electrodes.labels);
title("ICA Denoised NewData3");
xlim('tight');
grid minor
