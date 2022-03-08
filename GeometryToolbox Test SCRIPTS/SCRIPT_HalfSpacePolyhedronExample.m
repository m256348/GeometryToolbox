%% SCRIPT_HalfSpacePolyhedronExample
% This script defines the interior of a *convex* polyhedron using the
% intersection of negative half spaces.
%
%   NOTE: This script requires nCross from the TransformationToolbox
%
%   M. Kutzer, 25Feb2022, USNA

clear all
close all
clc

%% Create a random polygon
verts = rand(3,100);

% Make vertices convex
faces = convhull(verts(1,:),verts(2,:),verts(3,:));

%% Plot polygon
fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);
view(axs,3);

ptc = patch('Vertices',verts.','Faces',faces,'FaceColor','b',...
    'FaceAlpha',0.5,'EdgeColor','b','LineWidth',1);

%% Define half space for each segment of the polygon using Hessian Normal form
for k = 1:size(faces,1)
    % Define adjacent verticex pair
    p1 = verts(:,faces(k,1));
    p2 = verts(:,faces(k,2));
    p3 = verts(:,faces(k,3));

    % --- Define "oriented" plane using fitPlane --------------------------
    % NOTES:
    % (1) Creating a reliable line orientation requires nCross from the
    %     TransformationToolbox
    abc(k,:) = fitPlane([p1,p2,p3]);
    % ---------------------------------------------------------------------
end

%% Check half spaces for collisions
xlim(axs,[-1,2]);
ylim(axs,[-1,2]);
zlim(axs,[-1,2]);

% Define plot objects for interior and exterior points
plt_in  = plot3(axs,nan,nan,nan,'*r');
plt_out = plot3(axs,nan,nan,nan,'*g');

while true
    tic
    % Define test point
    p_tst = 3*rand(3,1) - 1;

    % Append a 1 to p
    p_tst(end+1) = 1;

    % Check point relative to line
    tst = abc * p_tst; 

    % Check if point is in each negative half space
    tf = tst <= 0; % See who's negative half space we are in
    
    % Check if point if in the intersection of negative half spaces
    if nnz(tf) ~= numel(tf)
        % Point *is not* in the intersection of negative half spaces
        x = get(plt_out,'XData');
        y = get(plt_out,'YData');
        z = get(plt_out,'ZData');
        x(end+1) = p_tst(1);
        y(end+1) = p_tst(2);
        z(end+1) = p_tst(3);
        set(plt_out,'XData',x,'YData',y,'ZData',z);
    else
        % Point *is* in the intersection of negative half spaces
        x = get(plt_in,'XData');
        y = get(plt_in,'YData');
        z = get(plt_in,'ZData');
        x(end+1) = p_tst(1);
        y(end+1) = p_tst(2);
        z(end+1) = p_tst(3);
        set(plt_in,'XData',x,'YData',y,'ZData',z);
    end
    toc
    drawnow
end

