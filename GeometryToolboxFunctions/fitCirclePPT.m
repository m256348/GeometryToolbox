function [cfit,Xint_cfit] = fitCirclePPT(X1,X2,abc,varargin)
% FITCIRCLEPPT fits a circle to two points and a line assuming the line is
% tangent and coincident to the circle.
%   [cfit,Xint_cfit] = fitCirclePPT(X1,C2,abc)
%   ___ = fitCirclePPT(___,debug)
%   ___ = fitCirclePPT(___,ZERO)
%
%   Input(s)
%       X1 - 2x1 array defining first 2D point on the circle
%       X2 - 2x1 array defining second 2D point on the circle
%    abc - 1x3 array defining the coefficients to a line that is both
%           coincident and tangent to the circle
%   debug - [OPTIONAL] logical scalar indicating whether or not to show 
%           debug plots. Default value is false.
%    ZERO - [OPTIONAL] positive value that is sufficiently close to zero 
%           to be assumed zero (e.g. ZERO = 1e-8). If ZERO is not 
%           specified or if ZERO = [], a default value is used.
%
%   Output(s)
%        cfit - 2-element structured array containing the following fields
%           cfit.Center - 3x1 center of the circle (z-coordinate is 0)
%           cfit.Normal - 3x1 normal to the circle ([0; 0; 1])
%           cfit.Radius - radius of the circle
%   Xint_cfit - 2-element cell array defining the points of intersection
%               between the lines and the circle
%       Xint_cfit{i}(:,j) - x/y intersection of circle i with line j
%
%   NOTE: cfit fields are to retain compatibility with fitCircle.m
%
%   See also fitCircle fitPlane
%
%   B.A.Todd, M. Kutzer, 25OCT2024, USNA

%% Check input(s)
narginchk(3,5);

if numel(X1) ~= 2 || ~isnumeric(X1)
    error('Point on circle must be 2D.');
end

if numel(X2) ~= 2 || ~isnumeric(X2)
    error('Point on circle must be 2D.');
end

if numel(abc) ~= 3 || ~isnumeric(abc)
    error('Line must be defined using three constants.');
end

%% Reshape input(s)
X1 = reshape(X1,[],1);
X2 = reshape(X2,[],1);
abc = reshape(abc,1,[]);

a = abc(1);
b = abc(2);
c = abc(3);

% Set defaults
debug = true;
ZERO = 1e-8;
if nargin > 3
    % Parse values
    for i = 1:numel(varargin)
        if numel(varargin{i}) > 1
            error('Optional values must be scalars.');
        end
        if numel(varargin{i}) == 0
            continue
        end

        if islogical(varargin{i})
            debug = varargin{i};
        else
            ZERO = varargin{i};
        end
    end
end

if ZERO < 0
    error('ZERO must be a positive value.');
end

%% Initialize output
cfit.Center = [];
cfit.Normal = [0;0;1];
cfit.Radius = [];

%% Check for special case
tfSpecialCase = false(1,2);

% Check if either point lies on line
% Point 1 check
if abc*[X1;1] < ZERO
    tfSpecialCase(1) = true;
end

% Point 2 check
if abc*[X2;1] < ZERO
    tfSpecialCase(2) = true;
end

% Point lies on both Line 1 and Line 2
if all(tfSpecialCase,'all')
    msg = sprintf([...
        'Tangent line intersects both points\n',...
        'The resultant circle has a radius of 0.']);
    warning(msg);
    cfit.Center = [X1; 0];
    cfit.Radius = 0;
    Xint_cfit{1} = X1;
    return
end

%% Create debug plot
if debug
    % Create figure and axes
    fig = figure('Name','fitCirclePTT.m, debug = true');
    axs = axes('Parent',fig,'DataAspectRatio',[1 1 1],'NextPlot','add');

    % Plot points
    plt_X1 = plot(axs,X1(1),X1(2),'ok','MarkerFaceColor','k','DisplayName', 'X1');
    plt_X2 = plot(axs,X2(1),X2(2),'ok','MarkerFaceColor','k', 'DisplayName', 'X2');
end


%% Find perpendicular to line at the intersection point
% HELP: Isn't this hte equation for tangent?
ab_p = nCross([a,b]);

if tfSpecialCase(1)
    c_p = -[ab_p(1), ab_p(2)]*X1;
end

if tfSpecialCase(2)
    c_p = -[ab_p(1), ab_p(2)]*X2;
end

abc_p = [ab_p(1), ab_p(2), c_p];
%% Find segment connecting two points
pts = [X1, X2];
seg = fitLine(pts); % seg = [a b c]

%% Find midpoint of segment
midseg = [(X1(1)+X2(1))/2; (X1(2)+X2(2))/2];

%% Find perpendicular line from midpoint on segment
% Calculate the segment vector
seg = X2 - X1;  % Vector from X1 to X2

if seg(1) == 0 && seg(2) ~= 0
    % Segment is vertical; perpendicular line is horizontal
    a_perp_seg = 0;  % Coefficient for x
    b_perp_seg = 1;  % Coefficient for y
    c_perp_seg = -midseg(2);  % y-intercept
elseif seg(2) == 0 && seg(1) ~= 0
    % Segment is horizontal; perpendicular line is vertical
    a_perp_seg = 1;  % Coefficient for x
    b_perp_seg = 0;  % Coefficient for y
    c_perp_seg = -midseg(1);  % x-intercept
else
    % Calculate the slope of the original segment
    slope = seg(2) / seg(1);  % Change in y over change in x
    % Calculate the slope of the perpendicular line
    perp_slope = -1 / slope;  % Negative reciprocal

    % Coefficients of the perpendicular line in the form ax + by + c = 0
    a_perp_seg = perp_slope;   % Coefficient for x
    b_perp_seg = -1;            % Coefficient for y (standard form)
    c_perp_seg = -(a_perp_seg * midseg(1) + b_perp_seg * midseg(2));  % Calculate c
end

abc_p_s = [a_perp_seg, b_perp_seg, c_perp_seg];
%% Calculate intersect between both perpendicular lines
% -> Point of intersection is the same regardless of line orientation
Xint = intersectLineLine(abc_p(1,:),abc_p_s(1,:));

% Define distance between points
dX = norm(X1-Xint);
xx = Xint(1) + dX*[-1,1];
yy = Xint(2) + dX*[-1,1];

%% Package cfit, Xint_cfit
cfit.Center = [Xint; 0];
cfit.Radius = dX;
cfit.Normal = [0; 0; 1];
Xint_cfit{1} = [X1, X2];

%% Update debug
if debug
    % Intersection point
    plt_Xint = plot(axs, Xint(1), Xint(2), 'or', 'MarkerFaceColor', 'r', 'DisplayName', 'Intersection');
    hold on; % Ensure all plots are on the same axes
    % Midpoint
    plt_mid = plot(axs, midseg(1), midseg(2), 'ob', 'MarkerFaceColor', 'b', 'DisplayName', 'Midpoint'); 
    % Segment between X1 and X2
    plt_midseg = plot(axs, [X1(1), X2(1)], [X1(2), X2(2)], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Segment'); 
    % Original tangent line
    plt_abc = line(xlim, (abc(1) * xlim + abc(3)) / -abc(2), 'Color', 'm', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Original Tangent Line'); 
    % Perpendicular line from X1
    plt_perpX1 = line(xlim, (abc_p(1) * xlim + abc_p(3)) / -abc_p(2), 'Color', 'g', 'LineStyle', '-.', 'LineWidth', 1.5, 'DisplayName', 'Perpendicular from X1'); 
    % Perpendicular line from midpoint
    plt_perpseg = line(xlim, (abc_p_s(1) * xlim + abc_p_s(3)) / -abc_p_s(2), 'Color', 'c', 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Perpendicular from Midpoint'); 
    % Circle plot (assuming plotCircle returns a plot handle)
    plt_circ = plotCircle(cfit); % Update to return a handle for the circle plot
    % Create the legend with correct handles and entries
    lgnd = legend('show'); % Automatically gathers all 'DisplayName' properties
    lgnd.Location = 'best'; % Automatically finds the best location for the legend
    hold off; % Release the hold on the axes
end
end
