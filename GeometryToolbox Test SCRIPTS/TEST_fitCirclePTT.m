%% TEST_fitCirclePTT
% Test fitCirclePTT
%
%   M. Kutzer, 14May2024, USNA
clear all
close all
clc

%% Create random data set
pnts = 10*(rand(2,5)-0.5);

X = pnts(:,1);
abc1 = fitPlane(pnts(:,2:3));
abc2 = fitPlane(pnts(:,4:5));
%load('TestCase_fitCirclePTT.mat');

%% Test function
[cfit,Xint_cfit] = fitCirclePTT(X,abc1,abc2,true);

