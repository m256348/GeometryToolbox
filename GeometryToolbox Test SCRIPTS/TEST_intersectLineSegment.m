%% TEST_intersectLineSegment
% Test intersectLineSegment
%
%   M. Kutzer, 19Sep2024, USNA

%% Create figure and axes
fig = figure;
axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);

%% Define text points
n = 2;
X1 = rand(n,2);
X2 = rand(n,2);

abc = fitLine(X1);
seg = fitSegment(X2);

%% 