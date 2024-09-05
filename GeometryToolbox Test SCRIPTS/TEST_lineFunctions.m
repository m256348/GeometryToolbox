%% TEST_lineFunctions
clear all
close all
clc

%% Generate random points
pnts = 10*(rand(2,2) - 0.5);
abc = fitPlane(pnts);

fig = figure('Name','TEST_lineFunctions');
axs = axes('Parent',fig,'DataAspectRatio',[1 1 1],'NextPlot','add');
plt = plot(axs,pnts(1,:),pnts(2,:),'xk');
pltLine = plotLine(axs,abc);