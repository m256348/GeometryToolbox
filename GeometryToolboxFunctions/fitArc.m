function [afit,X_sort] = fitArc(X)
% FITARC fits a circular arc to a three-dimensional set of data.
%   afit = FITCIRCLE(X) fits a *single* arc to set of N three-dimensional 
%   points.
%
%   [afit,X_sort] = FITCIRCLE(X)
%
%   Input(s)
%       X    - 3xN array containing points
%
%   Output(s)
%       afit - structured array containing the following fields
%           afit.Center    - 3x1 center of the arc
%           afit.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.Radius    - radius of the arc
%           afit.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection. 
%     X_sort - returns the array of X values sorted along the progression 
%              of the arc
%
%   See also plotArc plotArcs arcs2xyy
%
%   M. Kutzer, 23Sep2021, USNA

% TODO - add mean reprojection error to output

%% Check Inputs
narginchk(1,1);

% Check points
[m,n] = size(X);
if m ~= 3
    error('Specified points must be provided as a 3xN array.');
end

%% Fit circle
cfit = fitCircle(X);

%% Project points to circle
X_w = proj2circle(cfit,X);
X_w(4,:) = 1;

%% Define rigid body transform
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

%% Reference points to body-fixed frame
X_c = invSE(H_c2w)*X_w;

%% Define angles for all points
thetas_p2p = atan2(X_c(2,:),X_c(1,:)); % [-pi,pi]
thetas_z2p = wrapTo2Pi(thetas_p2p);    % [0,2*pi]

% Sort candidate arc bounds
[thetas_p2p,idx_p2p] = sort(thetas_p2p);
[thetas_z2p,idx_z2p] = sort(thetas_z2p);

% TODO - Check if this is a reasonable way to accomplish ordering 
if max(abs(diff(thetas_p2p))) < max(abs(diff(thetas_z2p)))
    thetas = thetas_p2p;
    idx = idx_p2p;
else
    thetas = thetas_z2p;
    idx = idx_z2p;
end

% rad2deg( thetas_p2p(end) - thetas_p2p(1) )
% rad2deg( thetas_z2p(end) - thetas_z2p(1) )

%% Package outputs
afit.Center = reshape(cfit.Center,[],1);
afit.Rotation = H_c2w(1:3,1:3);
afit.Radius = cfit.Radius;
afit.AngleLims = [thetas(1),thetas(end)];

X_sort = X(:,idx);