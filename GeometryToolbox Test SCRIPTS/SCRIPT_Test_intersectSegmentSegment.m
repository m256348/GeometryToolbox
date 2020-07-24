%% SCRIPT_Test_intersectSegmentSegment
% This script tests the segment/segment intersection code to see if/when it
% fails
%
%   M. Kutzer, 20Jul2020, USNA
clear all
close all
clc

%% Create figure for visualizing results
h.fig = figure('Name','Test intersectSegmentSegment');
h.axs = axes('Parent',h.fig);
hold(h.axs,'on');
daspect(h.axs,[1 1 1]);
colors = 'rb';
for i = 1:2
    h.seg(i) = plot(h.axs,nan,nan,colors(i),'LineWidth',1.5,'Marker','o');
    h.ext(i) = plot(h.axs,nan,nan,'k:','LineWidth',1);
end
set(h.seg(1),'MarkerSize',10);
h.pnt = plot(h.axs,nan,nan,'*k');

%% Test general case
for i = 1:20
    % Generate random end-points for two segments
    for j = 1:2
        % Generate random points
        p{j} = 100*rand*2*(rand(2,2) - 0.5);
    end
    
    checkANDdisp(p,h);
    pause;
    
    % Generate special case(s)
    dp1 = diff(p{1},1,2);
    theta0 = atan2(dp1(1),dp1(2));
    m = norm( diff(p{2},1,2) ); % Magnitude of second segment
    for theta = linspace(theta0,theta0+2*pi,6)
        p{2} = [p{1}(:,1), p{1}(:,1)+[m*sin(theta); m*cos(theta)]];
        checkANDdisp(p,h);
        pause
    end
    
    for dm = 2*m*(rand(1,5)-0.5)
        for theta = linspace(theta0,theta0+2*pi,6)
            p{2} = [p{1}(:,1), p{1}(:,1)+ m*[sin(theta); cos(theta)]] + dm*[sin(theta0); cos(theta0)];
            checkANDdisp(p,h);
            pause
        end
    end
end

%% Internal function(s)
function checkANDdisp(p,h)
% Plot segments
for j = 1:2
    set(h.seg(j),'XData',p{j}(1,:),'YData',p{j}(2,:));
end
% Check for intersections for random segments
[intEE,intEV,intVV,intPt] = intersectSegmentSegment(p{1}, p{2});
% Plot intersection (if applicable)
if ~isempty(intPt)
    set(h.pnt,'XData',intPt(1,:),'YData',intPt(2,:),'Visible','on');
else
    set(h.pnt,'Visible','off');
end
% Extend segment if intersection is off-segment
if ~any([intEE,intEV,intVV]) && ~isempty(intPt)
    for j = 1:2
        set(h.ext(j),'XData',[p{j}(1,1),intPt(1,1)],'YData',[p{j}(2,1),intPt(2,1)],'Visible','on');
    end
else
    set(h.ext,'Visible','off');
end
% Show result in title
msg = sprintf('Edge/Edge: %d | Edge/Vert: %d | Vert/Vert: %d',intEE,intEV,intVV);
title(h.axs,msg);
drawnow;
end