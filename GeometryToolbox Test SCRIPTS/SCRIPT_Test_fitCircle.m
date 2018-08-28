%% SCRIPT_Test_fitCircle
clear all
close all
clc

%% Define points
% Define radius
r = 10;
% Define center
C = [1; 2; 3];
% Define rotation
H = Tx(C(1))*Ty(C(2))*Tz(C(3))*Rx(3)*Ry(2)*Rx(5)*Rz(3);
% Define portion of circle and number of points
n = 50;
frac = 0.25;
theta = linspace(0,frac*(2*pi),n);

% Define points
Xin = [r*cos(theta); r*sin(theta)];
Xin(3,:) = 0;

noiseLevel = 0;%(1/100)*r;
Xnoise = Xin + noiseLevel*(rand(3,n)-0.5);
Xnoise(4,:) = 1;

X = H*Xnoise;
X(4,:) = [];

%% Fit circle
[cfit,meanError] = fitCircle(X)
