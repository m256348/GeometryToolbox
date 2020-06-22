function h = plotCircle(varargin)
% PLOTCIRCLE plots a circle defined in 3D.
%   h = PLOTCIRCLE(cfit) plots a circle with a default of 100 equally
%   spaced points.
%
%   h = PLOTCIRCLE(cfit,N) plots a circle with a user-specified number of
%   points.
%
%   h = PLOTCIRCLE(axs,___) allows the user to specify the parent of the
%   plot object.
%
%   Inputs:
%        axs - *optional* handle of the parent of the plotted circle
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%          N - *optional* number of points used to plot the circle
%
%   Outputs:
%       h - patch object handle for the plotted circle
%
%   M. Kutzer, 20Jun2020, USNA

%% Check inputs
narginchk(1,3);

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

% Define circle
if nargin >= nargin_i
	cfit = varargin{nargin_i};
    nargin_i = nargin_i+1;
end

% Define number of points
if nargin >= nargin_i
    N = varargin{nargin_i};
else
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

%% Patch result
%h = plot3(X_w(1,:),X_w(2,:),X_w(3,:),'Parent',axs,'LineWidth',2);
c.Vertices = X_w(1:3,:).';
c.Faces = 1:N;
h = patch(c,'FaceColor','None','EdgeColor','b','LineWidth',2);