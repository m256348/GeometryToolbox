%% SCRIPT_mvnpdf_ConfidenceInterval_ND
clear all
close all
clc

%% Parameters
N = 1;
kSamples = 1000000;

%% Define random statistics
% Define scale terms
a = 10;
b = 5;

% Generate a random 2x2 matrix
Sigma = a*rand(N,N);
% Multiply by its transpose
Sigma = Sigma.' * Sigma;
% Add a small multiple of the identity matrix to make it positive definite
epsilon = 0.1; % Small value to ensure positive definiteness
Sigma = Sigma + epsilon * eye(N);

% -> Generate random 1x3 mean
mu = b*rand(1,N);

%% Generate random points
statSample = mvnrnd(mu,Sigma,kSamples);

%% Evaluate confidence interval (value > 0 and < 1)
% -> 95% confidence interval (i.e. 95% of data is contained)
%    ---> Use 0.95
%
% References
% [1] Standard deviation and coverage,
%     https://en.wikipedia.org/wiki/Normal_distribution

% Capture data
ALL_prcConfInterval_Sample = [];
ALL_prcConfInterval        = [];

Xs_o = statSample.';
for prcConfInterval = linspace(0.01,0.99,1000)

    % Define n-dimensional confidence interval
    efit = mvnpdfConfInterval(mu,Sigma,prcConfInterval);
    
    % Reference samples to mean and covariance
    Xs_e = invSO(efit.Rotation)*(Xs_o - efit.Center);
    
    % Find all points inside the confidence interval
    tfIn = ...
        sum(( Xs_e./repmat(efit.PrincipalRadii,1,kSamples) ).^2, 1) <= 1; 

    prcConfInterval_Sample = nnz(tfIn)/numel(tfIn);

    fprintf('Confidence Interval: [USER DEFINED] %5.2f%% | %5.2f%% [OF SAMPLES]\n',...
        prcConfInterval*100,prcConfInterval_Sample*100);

    % Keep all of the info
    ALL_prcConfInterval(end+1)        = prcConfInterval;
    ALL_prcConfInterval_Sample(end+1) = prcConfInterval_Sample;
end