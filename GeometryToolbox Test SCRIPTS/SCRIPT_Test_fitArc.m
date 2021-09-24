%% SCRIPT_Test_fitArc
%
%   M. Kutzer, 24Sep2021, USNA

clear all
close all
clc

%% Create data set
r = 25;
theta0 = deg2rad(5+350);
theta1 = deg2rad(20+350);
thetas = linspace(theta0,theta1,5);
X_c(1,:) = r*cos(thetas);
X_c(2,:) = r*sin(thetas);
X_c(3,:) = 0;

%% Randomize data order
X_c = X_c(:,randperm(size(X_c,2)));

%% Add noise to data
err_c = 2*( rand(size(X_c)) - 0.5*ones(size(X_c)) );
%X_c = X_c + err_c;

%% Re-orient data
H_c2w = randSE(100);
X_c(4,:) = 1;
X_w = H_c2w*X_c;
X_w(4,:) = [];

%% Test fitArc
[afit,X_sort] = fitArc(X_w);

%% Plot result
fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);
view(axs,3);
axis(axs,'tight');

plt_X_w = plot3(axs,X_w(1,:),X_w(2,:),X_w(3,:),'xb');
plt_A = plotArc(axs,afit,500);
plt_X_s = plot3(axs,X_sort(1,:),X_sort(2,:),X_sort(3,:),'oc');

for i = 1:size(X_w,2)
    txt_w(i) = text(X_w(1,i),X_w(2,i),X_w(3,i),sprintf('%d_w',i),...
        'Parent',axs,'HorizontalAlignment','Left');
    txt_s(i) = text(X_sort(1,i),X_sort(2,i),X_sort(3,i),sprintf('%d_s',i),...
        'Parent',axs,'HorizontalAlignment','Right');
end

%% Test interpArc
X = interpArc(afit,50);
plt_i = plot3(axs,X(1,:),X(2,:),X(3,:),'.m','MarkerSize',12);