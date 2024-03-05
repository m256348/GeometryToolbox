function afit = fitArc2pntRad(X,r)
% FITARC2PNTRAD fits a circular arc to two, two dimensional points given a
% raduis.
%
%   afit = fitArc2pntRad(X,r)
%
%   Input(s)
%       X - 2x2 array containing points
%           pnt(:,1) - x/y coordinate of point 1
%           pnt(:,2) - x/y coordinate of point 2
%       r - scalar value defining radius
%
%   Output(s)
%       afit - structured array containing the following fields
%           afit.Center    - 3x1 center of the arc
%           afit.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.Radius    - radius of the arc
%           afit.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection.
%
%   NOTE:
%       This function returns 2D arcs meaning the center of the arc will
%       have a z-component of 0.
%
%   M. Kutzer, 04Mar2024, USNA

%% Check Input(s)
narginchk(2,2)

% Check points
[m,n] = size(X);
if m ~= 2 || n ~= 2
    error('Specified points must be provided as a 2x2 array.');
end

% Check radius
if numel(r) ~= 1

end