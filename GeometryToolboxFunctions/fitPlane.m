function varargout = fitPlane(pnts)
% FITPLANE fits the equation of an N-dimensional plane to a set of points.
%   pln = FITPLANE(pnts) fit the equation of an N-dimensional plane 
%   (c_1*x + c_2*y + c_3*z ... + c_{N+1) = 0) to a provided set of points. 
%   Coefficients are calculated in the Hessian normal form using singular 
%   value decomposition to determine the unit normal.
%
%       pnts - NxM array containing points where M >= N
%       pln - 1x(N+1) array containing coefficients for plane equation
%             General: [c_1,c_2,c_3,... c_{N+1}] such that
%                      c_1*x + c_2*y + c_3*z ... + c_{N+1) = 0
%           Line (2D): [a,b,c] such that a*x + b*y + c = 0
%          Plane (3D): [a,b,c,d] such that a*x + b*y + c*z + d = 0
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
%   23Sep2021 - Updated to include general N-dimensional planes

%% check inputs
narginchk(1,1);

[n,m] = size(pnts);
if n < 2
    error('Points must be specified in N-dimensions where N > 1.');
end
if m < n
    error('At least %d points must be specified to define a %d-dimensional plane.',n,n);
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

pln = zeros(1,numel(n)+1);
pln(1:numel(n)) = n;
pln(numel(n)+1) = p;

%% package outputs
if nargout > 0
    varargout{1} = pln;
end

if nargout > 1
    [~,meanError] = proj2plane(pln,pnts);
    varargout{2} = meanError;
end