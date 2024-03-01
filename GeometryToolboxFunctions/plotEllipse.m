function h = plotEllipse(varargin)
% PLOTELLIPSE plots an ellipse in 2D or 3D
%   h = PLOTELLIPSE(cfit) plots a ellipse with a default of 100 equally
%   spaced points.
%
%   h = PLOTELLIPSE(cfit,N) plots a ellipse with a user-specified number of
%   points.
%
%   h = PLOTELLIPSE(axs,___) allows the user to specify the parent of the
%   plot object.
%
%   Inputs:
%        axs - *optional* handle of the parent of the plotted ellipse
%       efit - structured array containing the following fields
%           [2D CASE]
%               efit.Center         - 2x1 center of the ellipse
%               efit.Rotation       - 2x2 rotation of the ellipse
%               efit.PrincipalRadii - 2x1 radii of each principal semi-axis
%           [3D CASE]
%               efit.Center         - 3x1 center of the ellipse
%               efit.Rotation       - 3x3 rotation of the ellipse
%               efit.PrincipalRadii - 2x1 radii of each principal semi-axis
%          N - *optional* number of points used to plot the ellipse
%
%   Outputs:
%       h - patch object handle for the plotted ellipse
%
%   See also fitEllipse interpEllipse
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
	efit = varargin{nargin_i};
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

%% Interpolate points on circle
X_w = interpEllipse(efit,N);

%% Patch result
%h = plot3(X_w(1,:),X_w(2,:),X_w(3,:),'Parent',axs,'LineWidth',2);
c.Vertices = X_w.';
c.Faces = 1:N;
h = patch(c,'FaceColor','None','EdgeColor','b','LineWidth',2);