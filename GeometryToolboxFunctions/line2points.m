function [xy,xyBnds] = line2points(abc,xx,yy,ZERO)
% LINE2POINTS defines two points on a line given x and y bounds.
%   xy = line2points(abc,xx,yy)
%
%   Input(s)
%      abc - 1x3 array defining the coefficients to a line that is both
%            coincident and tangent to the circle
%       xx - 1x2 array defining lower and upper x limits for points
%            returned
%       yy - 1x2 array defining lower and upper y limits for points
%            returned
%     ZERO - [OPTIONAL] scalar value that is very close to 0. Default value
%            is 1e-8.
%
%   Output(s)
%       xy - 2x2 array defining points on the line bounded by xx and yy
%           xy(:,1) - x/y coordinates of point 1
%           xy(:,2) - x/y coordinates of point 2
%
%   See also fitPlane plotLine intersectLineSegment intersectLineLine
%
%   M. Kutzer, 14May2024, USNA

% Updates
%   17Sep2024 - For three or more points, choose the farthest two
%   19Sep2024 - Added ZERO
%   19Sep2024 - Revised to use intersectLineSegment and intersectLineLine

debug = false;

%% Check input(s)
narginchk(3,4);

if numel(abc) ~= 3 || ~isnumeric(abc)
    error('Line must be defined using three constants.');
end

if numel(xx) ~= 2 || ~isnumeric(xx)
    error('Bounds for x must be defined using two values.');
end

if numel(yy) ~= 2 || ~isnumeric(yy)
    error('Bounds for y must be defined using two values.');
end

if nargin < 4
    % Set default value
    ZERO = 1e-8;
end

%% Setup debug plot
if debug
    fig = figure('Name','line2points, debug = true');
    axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
end

%% Ensure proper ordering
xx = sort(xx);
yy = sort(yy);

%% Define corner points
xyC = [...
    xx(1), xx(2), xx(2), xx(1);...
    yy(1), yy(1), yy(2), yy(2)];

% DEBUG PLOT
if debug
    plt_xyC = plot(axs,xyC(1,:),xyC(2,:),'oc');

    for i = 1:size(xyC,2)
        txt_xyC = text(axs,xyC(1,i),xyC(2,i),sprintf('%d',i),...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end
end

%% Define edge indices
edge_i = [...
    1,2;...
    2,3;...
    3,4;...
    4,1];

%% Define edge segments
xy = [];
xyBnds = [];
for i = 1:size(edge_i,1)
    % Fit segment
    seg = fitSegment(xyC(:,edge_i(i,:)));
    
    % Find the intersection between line and segment
    Xint = intersectLineSegment(abc,seg,ZERO);
    
    % SPECIAL CASE: Parallel, overlapping line and segment
    if size(Xint,2) == 2
        xy = Xint;
        xyBnds = Xint;
        break
    end

    % Only consider the single intersection case
    %   - Two intersections should be assounted for with adjacent edges
    if size(Xint,2) == 1
        xy = [xy, Xint];
    else
        abc2 = fitLine(xyC(:,edge_i(i,:)));
        Xint = intersectLineLine(abc,abc2);
        xyBnds = [xyBnds,Xint];
    end
    
    % DEBUG PLOT
    if debug
        plt_seg = plot(axs,seg(1,:)*[0,1;1,1],seg(2,:)*[0,1;1,1],'-c');
    end
end

% SPECIAL CASE: Line contains at least one segment end-point
m = size(xy,2);
if m > 2
    % Define distances between points
    d = inf(m,m);
    for i = 1:(m-1)
        for j = (i+1):m
            d(i,j) = norm( xy(:,i) - xy(:,j) );
        end
    end
    
    for i = 1:(m-2)
        [~,jj] = find(d == min(min(d)),1,'first');
        xyBnds = [xyBnds,xy(:,jj)];
        xy(:,jj) = [];
        d(jj,:) = [];
        d(:,jj) = [];
    end
end

% DEBUG PLOT
if debug
    % Plot xy
    plt_xy = plot(axs,xy(1,:),xy(2,:),'x'); 
    for j = 1:size(xy,2)
        txt_xy(j) = text(axs,xy(1,j),xy(2,j),sprintf('s_{%d}',j),...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    end

    % Plot xyBnds
    plt_xyBnds = plot(axs,xyBnds(1,:),xyBnds(2,:),'-'); 
    for j = 1:size(xy,2)
        txt_xy(j) = text(axs,xyBnds(1,j),xyBnds(2,j),sprintf('b_{%d}',j),...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    end
end
