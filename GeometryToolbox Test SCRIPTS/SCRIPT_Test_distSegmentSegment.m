%% SCRIPT_Test_distSegmentSegment
clear all
close all
clc

%% Create figure
fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
view(axs,3);
daspect(axs,[1 1 1]);
xlabel(axs,'x');
ylabel(axs,'y');
zlabel(axs,'z');

sg1 = plot(axs,0,0,'o-r');  % Segment 1
sg2 = plot(axs,0,0,'o-b');  % Segment 2
dst = plot(axs,0,0,'--k');  % Connection between segments
txt = text(axs,0,0,'');     % Connection distance

%% Test shortest distance between two segments
N = 2;
if N < 3
    view(axs,2);
end

ZERO = 1e-8;
while true
    mag = 1000 * rand(1,2);
    p = 2*rand(N,4) - ones(N,4);
    p1 = mag(1)*p(:,1);
    p2 = mag(1)*p(:,2);
    p3 = mag(2)*p(:,3);
    p4 = mag(2)*p(:,4);
    
    % Calculate distance between segments
    [d_mag,d,d1,d2] = distSegmentSegment(p1,p2,p3,p4);
    
    if N < 3
        p1(3) = 0;
        p2(3) = 0;
        p3(3) = 0;
        p4(3) = 0;
    end
    
    if ishandle(fig)
        set(sg1,...
            'XData',[p1(1),p2(1)],...
            'YData',[p1(2),p2(2)],...
            'ZData',[p1(3),p2(3)]);
        set(sg2,...
            'XData',[p3(1),p4(1)],...
            'YData',[p3(2),p4(2)],...
            'ZData',[p3(3),p4(3)]);
        set(dst,...
            'XData',[d1(1),d2(1)],...
            'YData',[d1(2),d2(2)],...
            'ZData',[d1(3),d2(3)]);
        set(txt,...
            'Position',mean([d1,d2],2),...
            'String',sprintf('%.2f',d));
    else
        break
    end
    
    % Check if discovered points are on each of the segments
    X1 = zeros(3,2);
    X2 = zeros(3,2);
    s = [0,1];
    for i = 1:3
        X1(i,:) = polyfit(s,[p1(i),p2(i)],1);
        X2(i,:) = polyfit(s,[p3(i),p4(i)],1);
    end
    s1 = pinv(X1(:,1))*(d1 - X1(:,2));
    s2 = pinv(X2(:,1))*(d2 - X2(:,2));
    
    % Check if point is on segment 1
    %if s1 < s(1) || s1 > s(2)
    if zeroFPError( s(1)-s1,ZERO ) > 0 || zeroFPError( s1-s(2),ZERO ) > 0
        fprintf('Closest point on Segment 1 is not on the segment, s1 = %f.\n',s1);
        break
    end
    % Check if point is on segment 2
    %if s2 < s(1) || s2 > s(2)
    if zeroFPError( s(1)-s2,ZERO ) > 0 || zeroFPError( s2-s(2),ZERO ) > 0
        fprintf('Closest point on Segment 2 is not on the segment, s2 = %f.\n',s2);
        break
    end
    % Check returned distance
    if zeroFPError( abs(norm(d1 - d2) - d_mag),ZERO ) > 0
        fprintf('Distance measure does not match contact point distance, |d1 - d2| - d = %f.\n',norm(d1 - d2) - d);
        break
    end
    
    drawnow
    pause(0.5);
    
    mag = mag + rand(1);
end