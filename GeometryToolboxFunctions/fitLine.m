function varargout = fitLine(pnts)
% FITLINE fits the equation of a 2-dimensional line to a set of points.
%   abc = FITLINE(pnts) fit the equation of a 2-dimensional line 
%   (c_1*x + c_2*y + c_{3) = 0) to a provided set of points. 
%   Coefficients are calculated in the Hessian normal form using singular 
%   value decomposition to determine the unit normal.
%
%   [..., meanError] = FITLINE(...) additionally returns the mean error 
%   of the distance of provided points and projections.
%
%   Input(s)
%       pnts - 2xM array containing points where M >= 2
%
%   Output(s)
%       abc - 1x3 array containing coefficients for line equation
%           Line (2D): [a,b,c] such that a*x + b*y + c = 0
%       meanError - mean error of the distance of provided points and their
%                   projections onto the line
%
%       NOTE: If the minimum number of points is provided, the line will
%             be "oriented" based on the order of points provided using the 
%             n-dimensional cross product (nCross) if the
%             TransformationToolbox is installed.
%
%   See also proj2line fitPlane
%
%   M. Kutzer, 05Sep2024, USNA

%% Check input(s)
narginchk(1,1);

[nD,mPnts] = size(pnts);
if nD ~= 2
    error('Points must be specified in 2-dimensions.');
end

%% Fit line
[varargout{1:nargout}] = fitPlane(pnts);