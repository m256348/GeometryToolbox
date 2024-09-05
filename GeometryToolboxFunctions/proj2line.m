function varargout = proj2line(abc,X)
% PROJ2LINE projects a set of points to a specified line
%   Xproj = PROJ2LINE(abc,X) projects the M points contained in the NxM
%   array X to the line specified as either a 1x(N+1) array or a 
%   planeModel object. 
%
%       abc - 1x3 array containing coefficients for line equation
%           Line (2D): [a,b,c] such that a*x + b*y + c = 0
%       X   - 2xM array containing points where M >= N
%
%   [..., meanError] = PROJ2LINE(...) additionally returns the mean 
%   error of the Euclidean distance between the provided points and 
%   their corresponding projections.
%
%   [..., Error] = PROJ2LINE(...) additionally returns errors for each
%   projected point.
%
%   See also fitLine, proj2plane
%
%   M. Kutzer, 05Sep2024, USNA

%% Check input(s)
narginchk(2,2);

[nD,mPnts] = size(X);
if nD ~= 2
    error('Points must be specified in 2-dimensions.');
end

if numel(abc) ~= 3
    error('Line must be defined as a 3-element array.');
end

%% Project to line
[varargout{1:nargout}] = proj2plane(abc,X);