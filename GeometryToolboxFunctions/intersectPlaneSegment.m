function pnt = intersectPlaneSegment(abcd,pnts)
% INTERSECTPLANSEGMENT finds the point of intersection, if it exists,
% between a plane and a segment.
%   pnt = INTERSECTPLANSEGMENT(abcd,pnts)
%
%   Inputs:
%       abcd - 1x4 array containing the coefficients of a plane
%       pnts - 3x2 array containing two 3D end-points of a segment
%
%   References:
%       [1] https://mathworld.wolfram.com/Line-PlaneIntersection.html
%
%   M. Kutzer, 12Jun2020, USNA

% Updates:
%   27May2021 - Finished function

% TODO - adjust function to eliminate warning, potentially add two outputs

%% Check inputs
narginchk(2,2);

if numel(abcd) ~= 4
    error('Plane must be defined using 4 coefficients.');
end
abcd = reshape(abcd,1,[]);

% TODO - check inputs

%% Fit parametric line where s \in [0,1] defines the segment
M = pnts*[0, 1; 1, 1]^(-1);

%% Solve for s
% M = [m1,m2];
m1 = M(:,1);
m2 = M(:,2);
% abcd = [abc,d];
abc = abcd(1:3);
d = abcd(4);

s = (-d - abc*m2)/(abc*m1);

if s >= 0 && s <= 1
    pnt = M*s;
else
    warning('Intersection is outside the segment.');
    pnt = M*s;
end