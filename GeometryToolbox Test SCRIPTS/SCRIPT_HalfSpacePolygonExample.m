%% SCRIPT_HalfSpacePolygonExample
% This script defines the interior of a polygon using an intersection of
% negative half spaces.
%
%   M. Kutzer, 24Feb2022, USNA

clear all
close all
clc

%% Create a random polygon
verts = rand(2,100);

% Make vertices convex
idx = convhull(verts(1,:),verts(2,:));
verts = verts(:,idx(1:end-1));

%% Plot polygon
fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);

faces = 1:size(verts,2);
ptc = patch('Vertices',verts.','Faces',faces,'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','LineWidth',1);

%% Define half space for each segment of the polygon using Hessian Normal form
idx = [faces,1]; % (Indices including wrap-around condition)
for i = 1:size(verts,2)
    % Define adjacent verticex pair
    p1 = verts(:,idx(i));
    p2 = verts(:,idx(i+1));
    
    % Define vector direction of vertex pair
    v = p2 - p1;
    v_hat = v./norm(v);
    
    % Define normal to vertex pair
    n_hat = [v_hat(2); -v_hat(1)];

    % Define offset
    p = -n_hat.'*p1;

    % Package line coefficients
    abc(i,:) = [n_hat.', p]; % Coefficients for our line ax + by + c = 0
    
    %{
    % fitPlane does not currently orient the line!
    % -> See TransformationToolbox nCross (currently doesn't support 
    %    single-element cross product for 2D vector) 
    abc(i,:) = fitPlane([p1,p2]);
    %}

    % Plot vertices and labels
    plot(axs,verts(1,idx(i)),verts(2,idx(i)),'*k');
    text(verts(1,idx(i)),verts(2,idx(i)),sprintf('p_{%d}',i));
end

%% Check half spaces for collisions
xlim(axs,[-1,2]);
ylim(axs,[-1,2]);

% Define plot objects for interior and exterior points
plt_in  = plot(axs,nan,nan,'*r');
plt_out = plot(axs,nan,nan,'*g');

while true
    tic
    % Define test point
    p_tst = 3*rand(2,1) - 1;

    % Append a 1 to p
    p_tst(3) = 1;

    % Check point relative to line
    tst = abc * p_tst; 

    % Check if point is in each negative half space
    tf = tst <= 0; % See who's negative half space we are in
    
    % Check if point if in the intersection of negative half spaces
    if nnz(tf) ~= numel(tf)
        % Point *is not* in the intersection of negative half spaces
        x = get(plt_out,'XData');
        y = get(plt_out,'YData');
        x(end+1) = p_tst(1);
        y(end+1) = p_tst(2);
        set(plt_out,'XData',x,'YData',y);
    else
        % Point *is* in the intersection of negative half spaces
        x = get(plt_in,'XData');
        y = get(plt_in,'YData');
        x(end+1) = p_tst(1);
        y(end+1) = p_tst(2);
        set(plt_in,'XData',x,'YData',y);
    end
    toc
    drawnow

    pause;
end

