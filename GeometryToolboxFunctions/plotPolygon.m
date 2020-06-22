function h = plotPolygon(varargin)
% PLOTPOLYGON creates a patch object representation of a polygon.
%   h = PLOTPOLYGON(Xv) plots a polygon with defined by vertices prescribed
%   in Xv.
%
%   h = PLOTPOLYGON(axs,___) allows the user to specify the parent of the
%   polygon.
%
%   Inputs:
%       axs - *optional* handle of the parent of the plotted polygon
%        Xv - 3xN array containing points that are the vertices of the
%        polygon.
%
%   Outputs:
%       h - patch object handle for the plotted polygon
%
%   M. Kutzer, 22Jun2020, USNA

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

% Define vertices
if nargin >= nargin_i
	Xv = varargin{nargin_i};
end

%% Plot polygon
p_struct.Vertices = Xv.';
p_struct.Faces = 1:size(Xv,2);

h = patch(p_struct,'Parent',axs,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','k');