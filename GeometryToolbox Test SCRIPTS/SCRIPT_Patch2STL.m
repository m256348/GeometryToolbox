%% SCRIPT_Patch2STL
% This script creates a patch object of a convex shape, uses the convex 
% hull to ensure that the result is a triangle based mesh, creates a
% triangulation, and writes and STL.
%
%   M. Kutzer, 23Sep2022, USNA

clear all
close all
clc

%% Create convex 3D patch
sfit = sphereModel([0,0,0,50]);

fig = figure;
axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
view(axs,3);
p_o = plotSphere(sfit);
set(p_o,'FaceColor','g','FaceAlpha',0.9,'EdgeColor','g');

%% Ensure patch is a triangle based mesh 
verts = p_o.Vertices;
%faces = p_o.Faces;     % Tetrahedral connectivity
k = convhull(verts(:,1),verts(:,2),verts(:,3));

%% Create triangulation and write STL
%TR = triangulation(faces, verts);
TR = triangulation(k, verts);

fname = 'test.stl';
stlwrite(TR,'test.stl','binary');

%% Read STL and plot result
[TR,fileformat] = stlread(fname);
verts = TR.Points;
faces = TR.ConnectivityList;

p_f = patch('Faces',faces,'Vertices',verts,'FaceColor','b','FaceAlpha',0.5);