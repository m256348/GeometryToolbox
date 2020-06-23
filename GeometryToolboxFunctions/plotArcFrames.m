function h = plotArcFrames(varargin)
% PLOTARCFRAMES plots frames defined using the normal and tangent to the
% arc.
%
% arc.
%   h = PLOTARCFRAMES(afit) plots a default of 5 frames defined using the
%   normal and tangent to the arc.
%
%   h = PLOTARC(afit,N) allows the user to specify the number of frames.
%
%   h = PLOTARCFRAME(axs,___) allows the user to specify the parent of the 
%   plot object(s).
%
%   Inputs:
%        axs - *optional* handle of the parent of the plotted circle
%       afit - structured array containing the following fields
%           afit.Center    - 3x1 center of the arc
%           afit.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.Radius    - radius of the arc
%           afit.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection.
%          N - *optional* number of points used to plot the arc
%
%   Outputs:
%       h - hgtransform object handle(s) used to visualize the coordinate
%       frames along the arc.
%
%   M. Kutzer, 23Jun2020, USNA

%% Check inputs
narginchk(1,2);

% Define parent
parentGiven = false;
nargin_i = 1;
if ishandle(varargin{nargin_i})
    switch get(varargin{nargin_i},'Type')
        case 'axes'
            parentGiven = true;
        case 'hgtransform'
            parentGiven = true;
    end
end

if parentGiven
    axs = varargin{nargin_i};
    nargin_i = nargin_i+1;
else
    axs = gca;
end

switch get(axs,'Type')
    case 'axes'
        hold(axs,'on');
end

% Define arc
if nargin >= nargin_i
	afit = varargin{nargin_i};
nargin_i = nargin_i+1;
end

% Define number of frames
if nargin >= nargin_i
    N = varargin{nargin_i};
else
    N = 5;
end

%% Check for empty afit
if isempty(afit)
    h = [];
    return
end

%% Define circle frame
H_c2w = eye(4);
H_c2w(1:3,1:3) = afit.Rotation;
H_c2w(1:3,4) = afit.Center;

%% Plot arc frames
R_c2w = H_c2w(1:3,1:3);
for i = 1:size(afit.AngleLims,1)
    theta = linspace(afit.AngleLims(i,1),afit.AngleLims(i,2),N);
    
    % Define z-direction
    z_hat = R_c2w(1:3,3);
    % Define x-directions (unit tangent)
    x_hats = R_c2w * [-sin(theta); cos(theta); zeros(1,N)];
    
    % Define body-fixed points
    X_c(1,:) = afit.Radius.*cos(theta);
    X_c(2,:) = afit.Radius.*sin(theta);
    X_c(3,:) = 0;
    X_c(4,:) = 1;
    % Define world-referenced points
    X_w = H_c2w*X_c;
    
    for j = 1:N
        x_hat = x_hats(:,j);
        y_hat = cross(z_hat,x_hat);
        Hj_c2w = eye(4);
        Hj_c2w(1:3,1:3) = [x_hat,y_hat,z_hat];
        Hj_c2w(1:3,4) = X_w(1:3,j);
        
        h(i,j) = triad('Parent',axs,'Matrix',Hj_c2w,'Scale',1,'LineWidth',2);
    end
end