function afit = intersectCirclePolygon(cfit,Xv)
% INTERSECTCIRCLEPOLYGON defines an arc representing the intersection
% between a planar polygon and a circle lying on the same plane.
%   afit = INTERSECTCIRCLEPOLYGON(cfit,Xv) 
%
%   Inputs:
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%         Xv - 3xN array containing vertices of the polygon
%
%   Outputs:
%       afit - structured array containing the following fields
%           afit.Center      - 3x1 center of the arc
%           afit.Rotation    - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.AngleLimits - Nx2 array containing the bounds of the 
%                              angles used to define the arc.
%                AngleLimits(i,:) - lower and upper bounds of the angle
%                                   defining the ith arc intersection. 
%
%   M. Kutzer, 22Jun2020, USNA

% TODO - make "zero" value something that the user can define
ZERO = 1e-6;

%% Check inputs
narginchk(2,2);

% TODO - actually check the inputs

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

hh = triad('Matrix',H_c2w);
plot(hh,X_c(1,:),X_c(2,:),'db');

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
    b = dot(A,B);
    c = dot(B,B) - cfit.Radius;
    
    if (b^2 - 4*a*c) < 0
        % Point is imaginary
    else
        % Intersect (1)
        s = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        if s >= 0 && s <= 1
            X_int = M*[s; 1];
            plot(hh,X_int(1,:),X_int(2,:),'*m');
            theta(end+1) = atan2(X_int(2),X_int(1));
        end
        % Intersect (2)
        s = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        if s >= 0 && s <= 1
            X_int = M*[s; 1];
            plot(hh,X_int(1,:),X_int(2,:),'*m');
            theta(end+1) = atan2(X_int(2),X_int(1));
        end
    end
end
