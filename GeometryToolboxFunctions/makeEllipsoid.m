function efit = makeEllipsoid(varargin)
% MAKEELLIPSOID returns a structured array representation of an ellipsoid
% given a variety of inputs.
%   efit = MAKEELLIPSOID(V,D) returns an ellipsoid centered at the origin
%   with a rotation defined by eigenvectors V and principal radii defined 
%   by eigenvalues contained in the diagonal matrix D.
%
%       V    - 3x3 array containing directions of principal semi-axes
%       D    - 3x3 array containing magnitudes of principal semi-axes along 
%           the diagonal *or* a 3x1 containing the magnitudes of principal
%           semi-axes.
%       efit - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%
%   efit = MAKEELLIPSOID(V,D,C) returns an ellipsoid centered at C
%   with a rotation defined by eigenvectors V and principal radii defined 
%   by eigenvalues contained in the diagonal matrix D.
%
%       C    - 3x1 array containing the center of the ellipsoid
%
%   efit = MAKEELLIPSOID(A) returns an ellipsoid defined by a 4x4 affine
%   transformantion.
%       NOTE: The affine transformation used to specify the ellipse can
%       only contain translations, rotations, and scaling. 
%
%   See also fitEllipsoid, proj2ellipsoid, inEllipsoid
%
%   M. Kutzer, 03Jan2018, USNA

%% Check Inputs
narginchk(1,3);

if nargin > 1
    V = varargin{1};
    D = varargin{2};
    [~,n] = size(D);
    if n == 1
        D = diag(D);
    end
    
    if nargin < 3
        C = zeros(3,1);
    else
        C = varargin{3};
    end
else
    A = varargin{1};
    sV = A(1:3,1:3);
    d = sqrt( sum(sV.^2,1) );
    V = sV./d;
    D = diag(d);
    C = A(1:3,4);
end

% Check dimensions of V, D, and C
[m,n] = size(V);
if ~(m == n && n == 3)
    error('V must be a 3x3 array.');
end

[m,n] = size(D);
if ~(m == n && n == 3)
    error('D must be a 3x3 or a 3x1 array.');
end

[m,n] = size(C);
if ~(m == 3 && n == 1)
    error('C must be a 3x1 array.');
end

% TODO - Check V for unit length vectors

%% Define parameters for ellipsoid
R = V;
R(:,3) = cross( V(:,1), V(:,2) );

s = sign( dot(R(:,3),V(:,3)) );
if s == 0
    error('Length of z-direction must be one.');
end
D(3,3) = s*D(3,3);

efit.Center = C;
efit.Rotation = R;
efit.PrincipalRadii = diag(D);