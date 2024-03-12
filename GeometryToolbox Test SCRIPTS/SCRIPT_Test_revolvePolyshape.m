%% SCRIPT_Test_revolvePolyshape
% Test revolvePolyshape.m
%
%   M. Kutzer, 12Mar2024, USNA
clear all
close all
clc

%% Define random polyshape
X = rand(2,20);
k = convhull(X(1,:),X(2,:));
k(end) = [];
X = X(:,k);

ps = polyshape(X(1,:),X(2,:));

%% Plot polyshape
fig = figure('Name','SCRIPT_Test_revolvePolyshape');
axs = axes('Parent',fig,'DataAspectRatio',[1 1 1],'NextPlot','add');
%view(axs,3);

plt = plot(axs,ps,'FaceColor','b','EdgeColor','b');

%% Define axis of revolution
rAxis = 2*rand(2,2) + 1;

plt_a = plot(axs,rAxis(:,1),rAxis(:,2),'r');
plt_a0 = plot(axs,rAxis(1,1),rAxis(1,2),'xr');

%% Revolve polyshape
n = 100;
ptcStruct = revolvePolyshape(ps,rAxis,n);

ptc = patch(ptcStruct,'Parent',axs,'FaceColor','g','EdgeColor','k',...
    'FaceAlpha',0.5);