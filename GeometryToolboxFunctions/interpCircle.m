function X = interpCircle(cfit,N)
% INTERPCIRCLE interpolates a 3D circle fit into a set number of equally
% spaced points.
%   X = INTERPCIRCLE(cfit) creates 100 equally spaced points around the
%   circle.
%
%   X = INTERCIRCLE(cfit,N) create N equally spaced points around the
%   circle.
%
%   Inputs:
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%          N - *optional* number of points used to plot the circle
%
%   Outputs:
%       X - 3xN equally spaced points around the circle defined in cfit
%
%   See also fitCircle plotCircle
%
%   M. Kutzer, 16Sep2021, USNA

%% Set default(s)
% Define default number of points
if nargin < 2
    N = 100;
end

%% Check parameters
% TODO - check parameters

%% Define circle frame
% z-direction
z_hat = reshape(cfit.Normal./norm(cfit.Normal),[],1);
% "other" direction
n_hat = z_hat([2,3,1]);
% x-direction
x_hat = cross(n_hat,z_hat);
x_hat = x_hat./norm(x_hat);
% y-direction
y_hat = cross(z_hat,x_hat);

H_c2w = eye(4);
H_c2w(1:3,1:3) = [x_hat,y_hat,z_hat];
H_c2w(1:3,4) = reshape(cfit.Center,[],1);

%% Define body-fixed points
theta = linspace(0,2*pi,N+1);
theta(end) = [];

X_c(1,:) = cfit.Radius.*cos(theta);
X_c(2,:) = cfit.Radius.*sin(theta);
X_c(3,:) = 0;
X_c(4,:) = 1;

%% Define world-referenced points
X_w = H_c2w*X_c;
X_w = X_w(1:3,:);

%% Package output
X = X_w;