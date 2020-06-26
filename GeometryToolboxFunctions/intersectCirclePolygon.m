function afit = intersectCirclePolygon(cfit,Xv,varargin)
% INTERSECTCIRCLEPOLYGON defines an arc representing the intersection
% between a planar polygon and a circle lying on the same plane.
%   afit = INTERSECTCIRCLEPOLYGON(cfit,Xv) 
%
%   afit = INTERSECTCIRCLEPOLYGON(cfit,Xv,id)
%
%   Inputs:
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%         Xv - 3xN array containing vertices of the polygon
%         id - 1x1 *optional* identifying number/flag this is set to an 
%              empty set if no value is specified.
%
%   Outputs:
%       afit - structured array containing the following fields
%           afit.Center    - 3x1 center of the arc
%           afit.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.Radius    - radius of the arc
%           afit.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%           afit.ID        - 1x1 identifying value (empty if no value is
%                            specified)
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection. 
%
%   M. Kutzer, 22Jun2020, USNA

% TODO - make "zero" value something that the user can define
ZERO = 1e-6;

debugON = false;

%% Check inputs
narginchk(2,3);

if nargin < 3
    id = [];
else
    id = varargin{1};
end

% TODO - actually check the inputs

%% Check for empty cfit
if isempty(cfit)
    afit = [];
    return;
end

%% Define the polygon plane
abcd = fitPlane(Xv);

%% Check to see if planes are the same
% Center of the circle must lie on the polygon face plane
Xc = reshape(cfit.Center,[],1);
Xc(4) = 1;
chk = abcd * Xc;

if abs(chk) > ZERO
    error('The circle and face must lie on the same plane.\n\ta*x_c + b*y_c + c*z_c + d = %.8f',chk);
end

% Unit normals must be the same 
nC_hat = reshape(cfit.Normal,[],1)./norm(cfit.Normal);
nV_hat = reshape(abcd(1:3),[],1)./norm(abcd(1:3));

chk = dot(nC_hat,nV_hat);
if abs( abs(chk) - 1 ) > ZERO
    error('The circle and face must lie on the same plane.\n\tdot("circle normal","polygon normal") = %.8f',chk);
end

%% Define circle frame
% z-direction
z_hat = reshape(cfit.Normal./norm(cfit.Normal),[],1);
% "other" direction
n_hat = z_hat([2,3,1]);
% x-direction
x_hat = cross(n_hat,z_hat);
x_hat = x_hat./norm(x_hat);
% y-direction
y_hat = cross(z_hat,x_hat);

H_c2w = eye(4);
H_c2w(1:3,1:3) = [x_hat,y_hat,z_hat];
H_c2w(1:3,4) = reshape(cfit.Center,[],1);

%% Transform polygon vertices to circle frame
X_w = Xv;
X_w(4,:) = 1;

X_c = invSE(H_c2w) * X_w;

if debugON
    % Plot inputs
    fig = figure('Name','intersectCirclePolygon.m');
    axs = axes('Parent',fig);
    hold(axs,'on');
    daspect(axs,[1 1 1]);
    plotCircle(axs,cfit);
    plotPolygon(axs,Xv);
    % Plot body-fixed coordinate frame
    hh = triad('Matrix',H_c2w);
    % Plot body-fixed vertices
    plot(hh,X_c(1,:),X_c(2,:),'db');
end

%% Fit segments to edges of polygon & find intersection(s)
theta = [];
n = size(X_c,2);
for i = 1:n
    j = i+1;
    if j > n
        j = 1;
    end
    
    M = fitSegment([X_c(1:2,i),X_c(1:2,j)]);

    % Let M = [A,B]
    A = M(:,1);
    B = M(:,2);
    
    % a*s^2 + b*s + c = 0
    a = dot(A,A);
    b = 2*dot(A,B);
    c = dot(B,B) - (cfit.Radius)^2;
    
    if (b^2 - 4*a*c) < 0
        % Ignore - Point is imaginary
    else
        % Intersect (1)
        s = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        if s >= 0 && s <= 1 % s \in [0,1]
            X_int = M*[s; 1];
            theta(end+1) = atan2(X_int(2),X_int(1));
            
            if debugON
                plot(hh,X_int(1,:),X_int(2,:),'*m');
            end
        end
        % Intersect (2)
        s = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        if s >= 0 && s <= 1 % s \in [0,1]
            X_int = M*[s; 1];
            theta(end+1) = atan2(X_int(2),X_int(1));
            
            if debugON
                plot(hh,X_int(1,:),X_int(2,:),'*m');
            end
        end
    end
end

%% Check if no intersection occurs
if isempty(theta)
    afit = [];
    return;
end

%% Keep arcs that are contained in the polygon
thetas = sort(theta);
%thetas(end+1) = wrapTo2Pi(thetas(1));  % ERROR!
thetas(end+1) = thetas(1) + 2*pi;       % Correction!

% TODO - speed this up be removing loop!
for i = 2:numel(thetas)
    % Define mid-point of arc for testing
    theta = linspace(thetas(i-1),thetas(i),3);
    xq(i-1) = cfit.Radius.*cos(theta(2));
    yq(i-1) = cfit.Radius.*sin(theta(2));
end

in = inpolygon(xq,yq,X_c(1,:),X_c(2,:));

idx = find(in);
AngleLims = [thetas(idx).', thetas(idx+1).'];

if debugON
    for i = 1:size(AngleLims,1)
        theta = linspace(AngleLims(i,1),AngleLims(i,2),100);
        xx = cfit.Radius.*cos(theta);
        yy = cfit.Radius.*sin(theta);
        plot(hh,xx,yy,'m','LineWidth',3);
    end
end

%% Package outputs
afit.Center = reshape(cfit.Center,[],1);
afit.Rotation = H_c2w(1:3,1:3);
afit.Radius = cfit.Radius;
afit.AngleLims = AngleLims;
afit.ID = id;