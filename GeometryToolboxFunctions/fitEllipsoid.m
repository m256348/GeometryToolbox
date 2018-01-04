function varargout = fitEllipsoid(X)
% FITELLIPSOID fits an ellispoid to a three-dimensional set of data.
%   efit = FITELLIPSOID(X) fits an ellipsoid to a set of N
%   three-dimensional points.
%   
%       X    - 3xN array containing points
%       efit - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%
%   [..., meanError] = FITELLIPSOID(...) additionally returns the mean
%   error of the distance of provided points and their corresponding 
%   projections to the fit.
%
%   See also proj2Ellipsoid, inEllipsoid
%
%   M. Kutzer, 03Jan2018, USNA

warning('THIS FUNCTION IS INCOMPLETE');

%% Check Inputs
narginchk(1,1);

% Check points
[m,n] = size(X);
if m ~= 3
    error('Specified points must be provided as a 3xN array.');
end
