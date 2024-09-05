%% TEST_fitCirclePTT
% Test fitCirclePTT
%
%   M. Kutzer, 14May2024, USNA
clear all
close all
clc

%% Create random data set
%load('TestCase_fitCirclePTT.mat');

% Generate random points
pnts = 10*(rand(2,5)-0.5);
% Select first point
X = pnts(:,1);
% Fit planes
abc1 = fitPlane(pnts(:,2:3));
abc2 = fitPlane(pnts(:,4:5));

% Define special cases
X1 = proj2line(abc1,X);
X2 = proj2line(abc2,X);
X12 = intersectLineLine(abc1,abc2);

%% Test function
% -> Free point
[cfit  ,Xint_cfit  ] = fitCirclePTT(X  ,abc1,abc2,true);
% -> Point on Line 1
[cfit1 ,Xint_cfit1 ] = fitCirclePTT(X1 ,abc1,abc2,true);
% -> Point on Line 2
[cfit2 ,Xint_cfit2 ] = fitCirclePTT(X2 ,abc1,abc2,true);
% -> Point on Line 1 & 2
[cfit12,Xint_cfit12] = fitCirclePTT(X12,abc1,abc2,true);

