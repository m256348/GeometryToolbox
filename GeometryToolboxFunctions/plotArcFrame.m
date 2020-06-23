function h = plotArcFrame(varargin)
% PLOTARCFRAME plots the body-fixed coordinate frame used to define the 3D
% arc.
%   h = PLOTARCFRAME(afit) plots the body-fixed coordinate frame used to
%   define the arc.
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
%
%   Outputs:
%       h - hgtransform object handle used to visualize the body-fixed
%       coordinate frame.
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

%% Plot arc frame
h = triad('Parent',axs,'Matrix',H_c2w,'Scale',(2/3)*afit.Radius,'LineWidth',2);