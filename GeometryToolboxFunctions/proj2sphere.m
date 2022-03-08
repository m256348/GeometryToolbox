function varargout = proj2sphere(sph,X)
% PROJ2SPHERE projects a set of points to the surface of a specified sphere
%   Xproj = PROJ2SPHERE(sph,X) projects the N points contained in the 3xN
%   array X to the surface of a sphere specified as either a 1x4 array or a
%   sphereModel object.
%   
%       X   - 3x4 array containing points 
%       sph - 1x4 array containing sphere parameters [a,b,c,r] such that 
%           (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2 *or* a valid sphereModel
%           object
%
%   [..., meanError] = PROJ2SPHERE(...) additionally returns the mean 
%   error of the Euclidean distance between the provided points and 
%   their corresponding projections.
%
%   [..., Error] = PROJ2SPHERE(...) additionally returns errors for each
%   projected point.
%
%   See also fitSphere, pcfitsphere, sphereModel
%
%   M. Kutzer, 21Dec2017, USNA 

% Updates
%   07Mar2022 - Corrected documentation issue(s)
%% Check Inputs
narginchk(2,2);

% Get sphere parameters
switch lower( class(sph) )
    case 'spheremodel'
        C = transpose( sph.Center );
        r = sph.Radius;
    case 'double'
        C = transpose( sph(1:3) );
        r = sph(4);
    case 'single'
        C = transpose( sph(1:3) );
        r = sph(4);
    case 'sym'
        C = transpose( sph(1:3) );
        r = sph(4);
    otherwise
        error('Sphere must be specified as a 1x4 array or a valid sphereModel object.');
end

% Check points
[M,~] = size(X);
if M ~= 3
    error('Specified points must be provided as a 3xN array.');
end

%% Define vectors
% -> Vector from the center of the sphere to the specified point
v_vec = X - C;
v_hat = v_vec./sqrt( sum( v_vec.^2, 1) );

% -> Vector from the center of the sphere to the projected point
r_vec = r.*v_hat;

% -> Vector from specified point to the projected point
d_vec = v_vec - r_vec;

%% Define projected point
Xproj = X - d_vec;

%% Package outputs
if nargout > 0
    varargout{1} = Xproj;
end

if nargout > 1
    D = sqrt( sum( d_vec.^2, 1) );
    meanError = mean(D);
    varargout{2} = meanError;
end

if nargout > 2
    varargout{3} = D;
end
