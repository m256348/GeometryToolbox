function varargout = fitPlane(pnts)
% FITPLANE fits the equation of a plane (a*x + b*y + c*z + d = 0) to a 
% set of points.
%   pln = FITPLANE(pnts) fit the equation of a plane 
%   (a*x + b*y + c*z + d = 0) to a provided set of points. Coefficients are 
%   calculated in the Hessian normal form using singular value 
%   decomposition to determine the unit normal.
%
%       pnts - 3xN array containing points
%       pln - 1x4 array containing coefficients for plane equation 
%           [a,b,c,d] such that a*x + b*y + c*z + d = 0
%
%   [..., meanError] = FITPLANE(...) additionally returns the mean error 
%   of the distance of provided points and projections.
%
%   References
%       [1] http://mathworld.wolfram.com/Plane.html
%       [2] http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
%
%   See also pcfitplane, proj2plane, inPlane
%
%   M. Kutzer 10July2015, USNA

% Updates
%   20Dec2017 - Updated to include mean error calculation

%% check inputs
narginchk(1,1);

[m,n] = size(pnts);
if m ~= 3
    error('Points must be specified in 3D');
end
if n < 3
    error('At least 3 points must be specified to define a plane');
end

%% eliminate/ignore non-finite values
[~,j] = find(~isfinite(pnts));
pnts(:,j) = [];

%% define plane
N = size(pnts,2);
cg = (1/N)*sum(pnts,2);

%TODO - eliminate use of svd twice by actually looking at the values in the
%second output of the first call
[u,~,~] = svd( bsxfun(@minus,pnts,cg) );
s = svd( bsxfun(@minus,pnts,cg) );
idx = find(s == min(s),1);

n = u(:,idx);   % unit normal
p = -n'*cg;     % intercept

pln = zeros(1,4);
pln(1:3) = n;
pln(4) = p;

%% package outputs
if nargout > 0
    varargout{1} = pln;
end

if nargout > 1
    [~,meanError] = proj2plane(pln,pnts);
    varargout{2} = meanError;
end