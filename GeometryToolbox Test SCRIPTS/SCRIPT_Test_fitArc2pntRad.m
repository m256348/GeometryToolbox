%% SCRIPT_Test_fitArc2pntRad
clear all
close all
clc

% Define random points
sc = 100;
X = sc*rand(2,2) - sc/2;

v = diff(X,1,2);
r = sign(rand(1,1) - 0.5)*( norm(v)/2 + rand(1,1)*norm(v,2) );

% Test function
afit = fitArc2pntRad(X,r)