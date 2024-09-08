clc
close all
clear all

%% Part 1
load('Ex1.mat');
fs = 200;
t = 0.0001 : 1/fs : 51.2;
figure('Name',"Part1")
for i = 1:3
    subplot(3,1,i);
    plot(t,EEG_Sig(i,:));
    if i == 1
        title("EEG Signal of Ex1");
    end
    xlabel('Time(s)');
    ylabel("Channel " + i);
    xlim('tight');
    grid minor
end

%% Part 2
figure('Name',"Part2")
scatter3(EEG_Sig(1,:),EEG_Sig(2,:),EEG_Sig(3,:));
title("Scatter3 of EEG Signal of Ex1");
xlabel('Channel1 axis');
ylabel('Channel2 axis');
zlabel('Channel3 axis');
grid minor
xlim('tight');
ylim('tight');
zlim('tight');

%% Part 3
clc
%%%%% PCA
zero_mean_EEG = EEG_Sig - mean(EEG_Sig,2); % Zero mean EEG
Cov_EEG = cov(zero_mean_EEG.'); % Covariance matrix of zero mean EEG
[V,Lambda] = eig(Cov_EEG); % D is eigen value matrix and V is eigen vectors

%%%%% Sorting components by energy
[~,perm] = sort(diag(Lambda),'descend'); 
V = V(:,perm);
Lambda = Lambda(perm,perm);

%%%%% Producing sources matrix
D = diag(diag(Lambda).^(-0.5))* V.';
Sources = D*zero_mean_EEG;

%%%%% Plottings
figure('Name',"Part3-1")
scatter3(EEG_Sig(1,:),EEG_Sig(2,:),EEG_Sig(3,:)); % Main signal
title("3D plot of EEG Signal and PCA vectors");
xlabel('Channel1 axis');
ylabel('Channel2 axis');
zlabel('Channel3 axis');
grid minor
xlim('tight');
ylim('tight');
zlim('tight');
plot3dv(V(:,1)); % First component
plot3dv(V(:,2)); % Second component
plot3dv(V(:,3)); % Third component

figure('Name',"Part3-2")
for i = 1:3
    subplot(3,1,i);
    plot(t,Sources(i,:));
    if i == 1
        title("Estimated Sources");
    end
    xlabel('Time(s)');
    ylabel("Component " + i);
    xlim('tight');
    grid minor
end

figure('Name',"Part3-3")
scatter3(Sources(1,:),Sources(2,:),Sources(3,:)); % PCA sources signal
title("3D plot of estimated sources");
xlabel('Component1 axis');
ylabel('Component2 axis');
zlabel('Component3 axis');
grid minor
xlim('tight');
ylim('tight');
zlim('tight');

%% Part 4
clc
close all
%%%%% Use the pca function
zero_mean_EEG = EEG_Sig - mean(EEG_Sig,2); % Zero mean EEG
[C,~,L] = pca(zero_mean_EEG.');


%%%%% Producing sources matrix
D = diag((L.').^(-0.5))* C.';
Sources_2 = D * zero_mean_EEG;

%%%%% Plottings
figure('Name',"Part4-1")
scatter3(EEG_Sig(1,:),EEG_Sig(2,:),EEG_Sig(3,:)); % Main signal
title("3D plot of EEG Signal and PCA vectors");
xlabel('Channel1 axis');
ylabel('Channel2 axis');
zlabel('Channel3 axis');
grid minor
xlim('tight');
ylim('tight');
zlim('tight');
plot3dv(C(:,1)); % First component
plot3dv(C(:,2)); % Second component
plot3dv(C(:,3)); % Third component

figure('Name',"Part4-2")
for i = 1:3
    subplot(3,1,i);
    plot(t,Sources_2(i,:));
    if i == 1
        title("Estimated Sources");
    end
    xlabel('Time(s)');
    ylabel("Component " + i);
    xlim('tight');
    grid minor
end

figure('Name',"Part4-3")
scatter3(Sources_2(1,:),Sources_2(2,:),Sources_2(3,:)); % PCA sources signal
title("3D plot of estimated sources");
xlabel('Component1 axis');
ylabel('Component2 axis');
zlabel('Component3 axis');
grid minor
xlim('tight');
ylim('tight');
zlim('tight');

%% Part 5
clc

%%%%% Performing SVD on zero mean signal
[U,S,V] = svd(zero_mean_EEG);

%%%%% Displaying relationship between U in SVD(USV.') and C in PCA
U
C

%%%%% Displaying relationship between S in SVD(USV.') and Lambda in PCA
disp("S^2 / (N-1)");
diag(S).^2 / (length(EEG_Sig(1,:)) - 1)
Lambda

