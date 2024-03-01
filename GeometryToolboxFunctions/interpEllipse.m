function X = interpEllipse(efit,N)
% INTERPELLIPSE interpolates a 2D or 3D ellipse fit into a set number of
% equally spaced points.
%   X = INTERPELLIPSE(efit) creates 100 equally spaced points around the
%   ellipse.
%
%   X = INTERPELLIPSE(efit,N) create N equally spaced points around the
%   ellipse.
%
%   Inputs:
%       efit - structured array containing the following fields
%           [2D CASE]
%               efit.Center         - 2x1 center of the ellipse
%               efit.Rotation       - 2x2 rotation of the ellipse
%               efit.PrincipalRadii - 2x1 radii of each principal semi-axis
%           [3D CASE]
%               efit.Center         - 3x1 center of the ellipse
%               efit.Rotation       - 3x3 rotation of the ellipse
%               efit.PrincipalRadii - 2x1 radii of each principal semi-axis
%          N - *optional* number of points used to plot the ellipse
%
%   Outputs:
%       X - 2xN or 3xN equally spaced points around the ellipse defined in 
%           efit
%
%   NOTE - Points are not equally spaced. Arc length parameterization is
%          unfinished.
%
%   See also fitEllipse plotEllipse
%
%   M. Kutzer, 01Mar2024, USNA

%% Check input(s)
narginchk(1,2);

if nargin < 2
    N = 100;
end

% Check efit
flds = {'Center','Rotation','PrincipalRadii'};
if ~all( isfield(efit,flds) )
    msg = 'Ellipse must be specified by a structured array containing the following fields: ';
    for i = 1:numel(flds)
        if i < numel(flds)
            msg = sprintf('%s%s, ',msg,flds{i});
        else
            msg = sprintf('%sand %s.',msg,flds{i});
        end
    end
    error(msg);
end

% Check dimensions
n = numel( efit.Center );
if ~( n == 2 || n == 3 )
    error('Center of ellipse must be defined as a 2x1 or 3x1.')
end

if ~ismatrix(efit.Rotation) || size(efit.Rotation,1) ~= n || size(efit.Rotation,2) ~= n
    error('For an elipse with an %dx1 center, the rotation must be %dx%d.',n,n,n);
end

if numel(efit.PrincipalRadii) ~= 2
    error('The principal radii of the ellipse must be specified as a 2x1.');
end

%% Define points
% Sort principal radii
[efit.PrincipalRadii,idx] = sort(efit.PrincipalRadii,'descend');

% TODO - arc length parameterization! 
% Define arc length of ellipse in Q1
a = efit.PrincipalRadii(1);
b = efit.PrincipalRadii(2);

%{
m = sqrt( 1-(b/a)^2 );
[~,E] = ellipke(m);
sq = a*E;
sf = 4*sq;

% Interpolate over arc length
s = linspace(0,sf,N+1);
s(end) = [];

% Initialize variables
tfQ = false(4,N);
sQ = cell(4,1);
phiQ = cell(4,1);
for i = 1:4
    tfQ(i,:) = (s >= (i-1)*sq) & (s < (i)*sq);
    sQ{i} = s(tfQ(i,:)) - (i-1)*sq;

    phiQ_i = atan( (a/b)*tan(sQ{i}) );

    disp( rad2deg( phiQ_i ) )

    tfPhi = phiQ_i < 0;
    phiQ_i(tfPhi) = phiQ_i(tfPhi) + pi/2;
    
    disp( rad2deg( phiQ_i ) )

    phiQ{i} = phiQ_i + (i-1)*(pi/2);
end
phi = [phiQ{1},phiQ{2},phiQ{3},phiQ{4}];
%}

phi = linspace(0,2*pi,N+1);
phi(end) = [];

% Define body-fixed coordinates
X_e(1,:) = a*cos(phi);
X_e(2,:) = b*sin(phi);

% Account for switched indices
if idx(1) > idx(2)
    X_e = flipud(X_e);
end

%% Rotate and translate points
if n > 2
    X_e(3,:) = 0;
end

R_e2o = efit.Rotation;
d_e2o = reshape(efit.Center,[],1);

X = R_e2o*X_e + d_e2o;