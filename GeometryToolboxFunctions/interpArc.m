function X = interpArc(afit,N)
% INTERPARC interpolates a 3D arc fit into a set number of equally
% spaced points.
%   X = INTERPARC(afit) creates 100 equally spaced points around the
%   arc.
%
%   X = INTERPARC(afit,N) create N equally spaced points around the
%   arc.
%
%   Inputs:
%       afit - structured array containing the following fields
%           afit.Center    - 3x1 center of the arc
%           afit.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.Radius    - radius of the arc
%           afit.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection. 
%          N - *optional* number of points used to plot the circle
%
%   Outputs:
%       X - 3xN equally spaced points around the arc defined in afit
%
%   See also fitArc plotArc
%
%   M. Kutzer, 24Sep2021, USNA

%% Set default(s)
% Define default number of points
if nargin < 2
    N = 100;
end

%% Check parameters
% TODO - check parameters

%% Define circle frame
H_c2w = eye(4);
H_c2w(1:3,1:3) = afit.Rotation;
H_c2w(1:3,4) = afit.Center;

X_w = [];
for i = 1:size(afit.AngleLims,1)
    theta = linspace(afit.AngleLims(i,1),afit.AngleLims(i,2),N);
    
    % Define body-fixed points
    X_c(1,:) = afit.Radius.*cos(theta);
    X_c(2,:) = afit.Radius.*sin(theta);
    X_c(3,:) = 0;
    X_c(4,:) = 1;
    
    % Define world-referenced points
    X_w = [X_w, H_c2w*X_c];
end
X = X_w(1:3,:);