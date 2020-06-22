function M = fitSegment(X)
% FITSEGMENT fits a segment to two N-dimensional points
%   M = fitSegment(X) creates a matrix of coefficients for a parametric
%   representation of a segment such that:
%       X(:,1) = M*[0; 1];  % s = 0
%       X(:,2) = M*[1; 1];  % s = 1
%
%   Inputs:
%       X - Nx2 array containing the N-dimensional end points of the 
%           segment
%
%   Outputs:
%       M - Nx2 array containing coefficients
%
%   M. Kutzer, 22Jun2020, USNA

%% Check inputs
narginchk(1,1);

if size(X,2) ~= 2
    error('Two points must be specified. X should be an Nx2 array.');
end

% TODO - check X to make sure that the points are unique

%% Fit coefficients
S = [0, 1; 1, 1];
M = X * (S^(-1));
