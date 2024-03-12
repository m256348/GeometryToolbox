function ptcStruct = revolvePolyshape(ps,rAxis,n)
% REVOLVEPOLYSHAPE revolves a polyshape around a designated axis
%
% ptcStruct = revolvePolyshape(ps,rAxis,n)
%
%   Input(s)
%       ps - polyshape or polyshapes (holes are ignored)
%    rAxis - 2x2 array specifying the axis of revolution
%           rAxis(1,:) - origin of the axis of revolution
%           rAxis(2,:) - second point specifying the direction of the axis
%        n - scalar defining number of points around the returned 
%
%   Output(s)
%       ptcStruct - structured array or arrays defining the vertices and 
%                   faces for a patch object
%
%   M. Kutzer, 12Mar2024, USNA

%% Check input(s)
% TODO - check inputs

%% Remove holes in polyshape(s)
for i = 1:numel(ps)
    ps(i) = rmholes(ps(i));
    X_o{i} = ps(i).Vertices.';
    X_o{i}(4,:) = 1;
end

%% Define frame for revolution
% NOTE: With respect to Frame o, the axis of revolution should be the
%       z-direction extending from the origin.
d_p2o = rAxis(1,:).';
d_p2o(3,:) = 1;

x_tilde = X_p{1}(1:2,1).' - d_p2o(1:2,:);
x_tilde(3,:) = 0;

y_hat = rAxis(2,:).' - d_p2o(1:2,:);
y_hat = y_hat./norm(y_hat);
y_hat(3,:) = 0;

z_hat = cross(x_tilde,y_hat);
z_hat = z_hat./norm(z_hat);

x_hat = cross(y_hat,z_hat);
x_hat = x_hat./norm(x_hat);

H_p2o = [x_hat,y_hat,z_hat,d_p2o; [0,0,0,1]];
H_o2p = invSE(H_p2o);

%% Revolve polyshapes
phi = linspace(0,2*pi,n+1);
phi(end) = [];
for i = 1:numel(X_p)
    % Reference polyshape to axis of revolution
    X_p = H_o2p*X_o{i};

    % Define number of vertices in single cross-section
    n = size(X_p,2);

    % Revolve cross section
    verts = [];
    faces = [];
    for j = 1:numel(phi)
        H_p2r = Ry(phi(j));
        X_r = H_p2r*X_p;
        verts = [verts,X_r];
        if j == numel(phi)
            faces = [faces,...
                [1:N  ] + (j-1)*N;...
                [2:N,1] + (j-1)*N;...
                [2:N,1] + (0  )*N;...
                [1:N  ] + (0  )*N];
        else
            faces = [faces,...
                [1:N  ] + (j-1)*N;...
                [2:N,1] + (j-1)*N;...
                [2:N,1] + (j  )*N;...
                [1:N  ] + (j  )*N];
        end
    end

    % Reference vertices to original frame
    verts = H_p2o*verts;

    % Package output
    ptcStruct(i).Vertices = verts(1:3,:).';
    ptcStruct(i).Faces = faces.';
end


