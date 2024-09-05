function [val,abcd] = evalCfit(cfit,X)
% EVALCFIT evaluates the function of a 3D circle defined by cfit.
%   val = evalCfit(cfit,X)
%
%   Input(s)
%        cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle (z-coordinate is 0)
%           cfit.Normal - 3x1 normal to the circle ([0; 0; 1])
%           cfit.Radius - radius of the circle
%           X - NxM array (N \in {2,3}) array defining points to evaluate.
%
%   Output(s)
%       val - Mx2 array defined as follows
%           val(:,1) - given Xp defined as the projection of X onto the plane
%                    of the circle
%               val(:,1) = sum( (Xp - cfit.Center).^2 ) - cfit.Radius^2
%           val(:,2) - given abcd defined as the plane containing the circle
%                    abcd = [cfit.Normal.', cfit.Normal.'*cfit.Center]
%               val(:,2) = abcd*[X; 1]
%       pln - 1x4 array containing coefficients for plane equation
%          Plane (3D): [a,b,c,d] such that a*x + b*y + c*z + d = 0
%
%   Summary of values contained in val(i,:). Note that "p" is used to
%   represent a positive value and "n" represents a negative value.
%       [0,0] - X(:,i) lies on the circle and on the plane containing the
%               circle
%       [0,p] - X(:,i) lies in the positive half space of the plane, and
%               projection of point on plane lies on the circle
%       ...
%
%   M. Kutzer, 05Sep2024, USNA

%% Check input(s)
narginchk(2,2);

% TODO - check cfit

[n,m] = size(X);
if n < 2 || n > 3
    error('Point(s) to evaluate must be either 2xM or 3xM.');
end

if n == 2
    % Assume z = 0 for 2D points
    X(3,:) = 0;
end

%% Define circle plane
abcd = [cfit.Normal.', cfit.Normal.'*cfit.Center];

%% Project points to plane
Xp = proj2plane(abcd,X);
val(:,2) = ( abcd*[X; ones(1,size(X,2))] ).';

%% Define values
val(:,1) = ( sum( (Xp - repmat(cfit.Center,1,size(Xp,2)) ).^2 , 1 ) ...
    - cfit.Radius^2 ).';