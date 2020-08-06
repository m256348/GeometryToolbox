function X = arcs2xyz(afits,ds)
% ARCS2XYZ converts one or more of arcs into a series of world-referenced
% x/y/z coordinates equally spaced using a fixed arc length between points.
%   X = ARCS2XYZ(afits) defines x/y/z coordinates using a default arc
%   length of 0.1.
%
%   X = ARCS2XYZ(afits,ds) defines x/y/z coordinates using a specified arc
%   length of ds. 
%
%   NOTE: This function assumes the arcs specified are continuous in the
%         order given.
%
%   Inputs:
%       afits - n-element structured array containing the following fields
%           afits(i).Center    - 3x1 center of the arc
%           afits(i).Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afits(i).Radius    - radius of the arc
%           afits(i).AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(j,:) - lower and upper bounds of the angle
%                                 defining the jth arc intersection. 
%          ds - *optional* change in arc length between discrete points
%
%   Outputs:
%       X - 3xM element array of x/y/z coordinates
%
%   See also uniqueArcs, arcAdjacency, arcAdjacencyCycles, plotArcs
%
%   M. Kutzer, 05Aug2020, USNA

%% Check inputs
narginchk(1,2);

if nargin < 2
    ds = 0.1;
end

%% Check for empty afit
if isempty(afits)
    X = [];
    return
end

%% Define world points (assuming piecewise arcs are continuous) 
X_w = [];
for k = 1:numel(afits)
    afit = afits(k);
    
    % Define circle frame
    H_c2w = eye(4);
    H_c2w(1:3,1:3) = afit.Rotation;
    H_c2w(1:3,4) = afit.Center;
    
    % Define change in theta
    %   s = \theta r
    %   ds = d\theta r
    %   d\theta = ds / r
    dtheta = ds/afit.Radius;
    
    for i = 1:size(afit.AngleLims,1)
        % Define theta
        theta = afit.AngleLims(i,1):dtheta:afit.AngleLims(i,2);
        
        % Check for problem child
        if numel(theta) < 1
            warning('No angles produced to visualize arc!');
            fprintf(' -> afits(%d).AngleLims(%d,:) = [%.8f, %.8f]\n',k,i,afit.AngleLims(i,:));
            fprintf(' -> ds = %.8f | dtheta = %.8f\n',ds,dtheta);
            fprintf(' -> Using angle limits...\n');
            theta = afit.AngleLims(i,:);
        end
        
        % Define body-fixed points
        X_c = [];   % Re-initialize X_c
        X_c(1,:) = afit.Radius.*cos(theta);
        X_c(2,:) = afit.Radius.*sin(theta);
        X_c(3,:) = 0;
        X_c(4,:) = 1;
        
        % Append world-referenced points
        X_w = [X_w, H_c2w*X_c];
    end
end

%% Format output
X = X_w(1:3,:);