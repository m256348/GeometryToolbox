%% SCRIPT_Patch2STL

sfit = sphereModel([0,0,0,50]);
p_o = plotSphere(sfit);

verts = p_o.Vertices;
%faces = p_o.Faces;     % Tetrahedral connectivity
k = convhull(verts(:,1),verts(:,2),verts(:,3));

%TR = triangulation(faces, verts);
TR = triangulation(k, verts);

fname = 'test.stl';
stlwrite(TR,'test.stl','binary');

%% Read STL and plot
[TR,fileformat] = stlread(fname);
verts = TR.Points;
faces = TR.ConnectivityList;

figure;
patch('Faces',faces,'Vertices',verts,'FaceColor','b','FaceAlpha',0.5);