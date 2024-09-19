function [xy,xyBnds] = line2points(abc,xx,yy)
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
%
%   Output(s)
%       xy - 2x2 array defining points on the line bounded by xx and yy
%           xy(:,1) - x/y coordinates of point 1
%           xy(:,2) - x/y coordinates of point 2
%
%   See also fitPlane plotLine
%
%   M. Kutzer, 14May2024, USNA

% Updates
%   17Sep2024 - For three or more points, choose the farthest two
%   19Sep2024 - Added ZERO

% TODO - make this an input & set improved default values
ZERO = 1e-8;

%% Check input(s)
narginchk(3,3);

if numel(abc) ~= 3 || ~isnumeric(abc)
    error('Line must be defined using three constants.');
end

if numel(xx) ~= 2 || ~isnumeric(xx)
    error('Bounds for x must be defined using two values.');
end

if numel(yy) ~= 2 || ~isnumeric(yy)
    error('Bounds for y must be defined using two values.');
end

%% Define points
% Define x = (-b/a)*y + (-c/a)
x_abc = @(y,abc)(-abc(2)/abc(1))*y + (-abc(3)/abc(1));

% Define y = (-a/b)*x + (-c/b)
y_abc = @(x,abc)(-abc(1)/abc(2))*x + (-abc(3)/abc(2));

% Define bounded segments by finding where line intersects with axes
% limits
% -> x-lower
xy(1,1) = xx(1);
xy(2,1) = y_abc(xx(1),abc);
% -> x-upper
xy(1,2) = xx(2);
xy(2,2) = y_abc(xx(2),abc);
% -> y-lower
xy(1,3) = x_abc(yy(1),abc);
xy(2,3) = yy(1);
% -> y-upper
xy(1,4) = x_abc(yy(2),abc);
xy(2,4) = yy(2);

%% Define unique points 
xy = unique(xy.','rows').';

%% Define "bounds"
% This term is loosly defined! 
% TODO - do this properly
xyBnds = [xy(:,1),xy(:,end)];

%% Find the points inside of the limits
tf = ...
    xy(1,:) >= xx(1)*(1-ZERO) & xy(1,:) <= xx(2)*(1+ZERO) & ...
    xy(2,:) >= yy(1)*(1-ZERO) & xy(2,:) <= yy(2)*(1+ZERO);

switch nnz(tf)
    case 1
        % Only one point is within the specified limits
        xy = [];
    case 2
        % Two points returned (expected result)
        xy = xy(:,tf);
    otherwise
        %{
        % Multiple points returned (keep first two)
        xy = xy(:,1:2);
        %}
        idxPairs = nchoosek(find(tf),2);
        for i = 1:size(idxPairs,1)
            dxy(i) = norm(xy(:,idxPairs(i,1)) - xy(:,idxPairs(i,2)));
        end
        idxPair = find(dxy == max(dxy),1,'first');
        xy = xy(:,idxPairs(idxPair,:));
end

