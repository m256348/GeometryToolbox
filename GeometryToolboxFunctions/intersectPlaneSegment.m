function [pnt,tfEndPoint] = intersectPlaneSegment(abcd,pnts,ZERO)
% INTERSECTPLANSEGMENT finds the point of intersection, if it exists,
% between a plane and a segment.
%   [pnt,tfEndPoint] = intersectPlaneSegment(abcd,pnts)
%   ___ = intersectPlaneSegment(abcd,pnts,ZERO)
%   
%   Inputs:
%       abcd - 1x4 array containing the coefficients of a plane
%       pnts - 3x2 array containing two 3D end-points of a segment
%       ZERO - [OPTIONAL] positive scalar value close to zero. Default
%              value is 1e-8.
%
%   Output(s)
%              pnt - 3xN array describing point(s) of intersection
%                  -> 3x0 (empty set) if no intersection exists
%                  -> 3x1 if a single point of intersection exists
%                  -> 3x2 if segment is contained in the plane
%       tfEndPoint - 1x2 logical array describing whether one of both end
%                    point(s) of the segment lie on the plane. 
%
%   References:
%       [1] https://mathworld.wolfram.com/Line-PlaneIntersection.html
%
%   M. Kutzer, 12Jun2020, USNA

% Updates:
%   27May2021 - Updated function
%   10Oct2024 - Completed function

%% Check inputs
narginchk(2,3);

if nargin < 3
    ZERO = 1e-8;
end

if numel(abcd) ~= 4
    error('Plane must be defined using 4 coefficients.');
end
abcd = reshape(abcd,1,[]);

if size(pnts,1) ~= 3 || size(pnts,2) ~= 2
    error('Points in 3D must be defined as a 3x2 array');
end

%% Set default outputs
pnt        = zeros(3,0);
tfEndPoint = false(1,2);

%% Fit parametric line where s \in [0,1] defines the segment
M = pnts*[0, 1; 1, 1]^(-1);

%% Check for line parallel with plane
if abs( dot( M(:,1), abcd(1:3).' ) ) < ZERO
    % Segment and plane are parallel
    if any( abs( abcd*[pnts; 1,1] ) < ZERO )
        pnt = pnts;
    end

    % Define end-point condition
    tfEndPoint = true(1,2);

    return
end

%% Solve for s
% M = [m1,m2];
m1 = M(:,1);
m2 = M(:,2);
% abcd = [abc,d];
abc = abcd(1:3);
d = abcd(4);

s = (-d - abc*m2)/(abc*m1);

% Check specific cases
% -> Point 1 on plane
if abs(s) <= ZERO
    pnt = pnts(:,1);
    tfEndPoint = [true,false];
    return
end
% -> Point 2 on plane
if abs(s-1) <= ZERO
    pnt = pnts(:,2);
    tfEndPoint = [false,true];
    return
end
% -> Neither point on plane
%if s >= 0 && s <= 1
if (s+ZERO) >= 0 && (s-ZERO) <= 1
    pnt = M*[s; 1];
    tfEndPoint = false(1,2);
end