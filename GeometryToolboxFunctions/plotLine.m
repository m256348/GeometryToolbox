function [plt,xy] = plotLine(varargin)
% PLOTLINE plots a line given a set of coefficients defining the line
%
%   [plt,xy] = plotLine(abc) plots a line bounded by axes limits
%
%   [plt,xy] = plotLine(abc,xx,yy) plots a line bounded by xx and yy
%
%   ___ = plotLine(axs,___) specifies the axes object
%
%   ___ = plotLine(___,NAME,VALUE) specifies name/value pairs to set
%   properties for the line
%
%   Input(s)
%      axs - [OPTIONAL] axes or parent object of the line to be plotted
%      abc - 1x3 array defining the coefficients to a line that is both
%            coincident and tangent to the circle
%       xx - [OPTIONAL] 1x2 array defining lower and upper x limits for
%            points returned
%       yy - [OPTIONAL] 1x2 array defining lower and upper y limits for
%            points returned
%     NAME - [OPTIONAL] character array specifying line object property
%    VALUE - [OPTIONAL] line object property setting corresponding to
%            the specified name
%
%   Output(s)
%       plt - line object for plotted line
%        xy - 2x2 array defining points on the line bounded by xx and yy
%           xy(:,1) - x/y coordinates of point 1
%           xy(:,2) - x/y coordinates of point 2
%
%   See also fitPlane line2points

%% Parse input(s)
if nargin < 1
    error('The line must be defined using a 3-element array.')
end

% Initialize values
axs = [];
abc = [];
xx = [];
yy = [];

% Parse variable input argument
i = 0;
while true
    i = i+1;
    if i > numel(varargin)
        break
    end

    % Parse axes/parent
    if numel(varargin{i}) == 1 && ishandle(varargin{i})
        if ~isempty(axs)
            warning('Multiple parent objects defined, using first value.');
            varargin(i) = [];
            i = i-1;
            continue;
        end
        axs = varargin{i};
        varargin(i) = [];
        i = i-1;
        % TODO - confirm that the object is a valid parent for a line
        continue
    end

    % Parse abc, xx, or yy
    if isnumeric(varargin{i})
        switch numel(varargin{i})
            case 3
                if ~isempty(abc)
                    warning('Multiple lines defined, using first value.');
                    varargin(i) = [];
                    i = i-1;
                    continue
                end
                abc = varargin{i};
                varargin(i) = [];
                i = i-1;
                continue
            case 2
                if isempty(xx)
                    xx = varargin{i};
                    varargin(i) = [];
                    i = i-1;
                    continue
                end
                if isempty(yy)
                    yy = varargin{i};
                    varargin(i) = [];
                    i = i-1;
                    continue
                end
                warning('Multiple x/y limites defined, using first value.');
                varargin(i) = [];
                i = i-1;
                continue
        end
    end

    % Check if a NAME is defined
    if ischar(varargin{i}) || isstring(varargin{i})
        break
    end
end

% Check if line is undefined
if isempty(abc)
    error('Line must be specified as a 3-element array');
end

% Set default values if undefined
if isempty(axs)
    axs = gca;
end

if isempty(xx)
    xx = xlim(axs);
end

if isempty(yy)
    yy = ylim(axs);
end

xx = reshape(xx,1,[]);
yy = reshape(yy,1,[]);
abc = reshape(abc,1,[]);

%% Define points
[xy,xyBnds] = line2points(abc,xx,yy);
if size(xy,2) < 2
    str = sprintf([...
        'The specified line is not visible in the axes limits:\n',...
        '\txlim = [%.4f,%.4f]\n',...
        '\tylim = [%.4f,%.4f]\n\n',...
        'Adjusting limits to:\n',...
        '\txlim = [%.4f,%.4f]\n',...
        '\tylim = [%.4f,%.4f]'],...
        xx,yy,...
        [min([xx,xyBnds(1,:)]),max([xx,xyBnds(1,:)])],...
        [min([yy,xyBnds(2,:)]),max([yy,xyBnds(2,:)])]);
    warning(str);

    xy = xyBnds;
end

%% Plot line
if isempty(varargin)
    plt = plot(axs,xy(1,:),xy(2,:),'-');
else
    plt = plot(axs,xy(1,:),xy(2,:),varargin{:});
end