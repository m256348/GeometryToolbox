%% SCRIPT_Test_intersectPlanePlane

fig = figure;
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);

colors = 'rb';
for i = 1:2
    x{i} = rand(3,3);
    plt(i) = plot3(x{i}(1),x{i}(2),x{i}(3),['.',colors(i)]);
	ptc(i) = patch('Vertices',x{i}.','Faces',[1:3],'FaceColor',colors(i),...
        'FaceAlpha',0.5);

    abcd{i} = fitPlane(x{i});
end

M = intersectPlanePlane(abcd{1},abcd{2});

S = [-2,2; 1 1];
xx = M*S;

lin = plot3(xx(1,:),xx(2,:),xx(3,:),'-k');


