function h = plotArcs(varargin)
% PLOTARCS plots multiple arcs with poinst equally spaced using constant 
% arc length.
%   h = PLOTARCS(afits) plots one or more arcs with a default change in 
%   arc length of 0.1.
%
%   h = PLOTARCS(afits,ds) plots one or more arcs with a user specified
%   change in arc length.
%
%   h = PLOTARCS(axs,___) allows the user to specify the parent of the 
%   plot object(s).
%
%   NOTE: This function assumes the arcs specified are continues in the
%         order given.
%
%   Inputs:
%         axs - *optional* handle of the parent of the plotted arc
%       afits - structured array containing the following fields
%           afits.Center    - 3x1 center of the arc
%           afits.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afits.Radius    - radius of the arc
%           afits.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection. 
%          ds - *optional* change in arc length between discrete points
%
%   Outputs:
%       h - N-element array of line objects.
%
%   M. Kutzer, 26Jun2020, USNA

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

% Define arc
if nargin >= nargin_i
	afits = varargin{nargin_i};
    nargin_i = nargin_i+1;
end

% Define change in arc length
if nargin >= nargin_i
    ds = varargin{nargin_i};
else
    ds = 0.1;
end

%% Check for empty afit
if isempty(afits)
    h = [];
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

%% Plot arc
% [0.224,1.000,0.078] - Neon Green
h = plot3(X_w(1,:),X_w(2,:),X_w(3,:),'Parent',axs,'Color',[0.224,1.000,0.078],'LineWidth',2);
