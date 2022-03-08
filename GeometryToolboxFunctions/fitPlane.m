function varargout = fitPlane(pnts)
% FITPLANE fits the equation of an N-dimensional plane to a set of points.
%   pln = FITPLANE(pnts) fit the equation of an N-dimensional plane 
%   (c_1*x + c_2*y + c_3*z ... + c_{N+1) = 0) to a provided set of points. 
%   Coefficients are calculated in the Hessian normal form using singular 
%   value decomposition to determine the unit normal.
%
%   [..., meanError] = FITPLANE(...) additionally returns the mean error 
%   of the distance of provided points and projections.
%
%   Input(s)
%       pnts - NxM array containing points where M >= N
%
%   Output(s)
%       pln - 1x(N+1) array containing coefficients for plane equation
%             General: [c_1,c_2,c_3,... c_{N+1}] such that
%                      c_1*x + c_2*y + c_3*z ... + c_{N+1) = 0
%           Line (2D): [a,b,c] such that a*x + b*y + c = 0
%          Plane (3D): [a,b,c,d] such that a*x + b*y + c*z + d = 0
%       meanError - mean error of the distance of provided points and their
%                   projections onto the plane
%
%       NOTE: If the minimum number of points is provided, the plane will
%             be "oriented" based on the order of points provided using the 
%             n-dimensional cross product (nCross) if the
%             TransformationToolbox is installed.
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
%   24Feb2022 - Added plane orientation using nCross

%% Check inputs
narginchk(1,1);

[nD,mPnts] = size(pnts);
if nD < 2
    error('Points must be specified in N-dimensions where N > 1.');
end

%% Eliminate/ignore non-finite values
[~,j] = find(~isfinite(pnts));
pnts(:,j) = [];

%% Check supplied points
if mPnts < nD
    error('At least %d finite points must be specified to define a %d-dimensional plane.',nD,nD);
end

%% Check for nCross is minimum set of points is provided
useCross = false;
if nD == mPnts
    if exist('nCross','file') == 2
        useCross = true;
    else
        warning('Unable to orient plane. Download and install the kutzer\TransformationToolbox for this functionality.')
    end
end

%% Define plane
if ~useCross
    % Fit nD plane to mPnts > nD OR if nCross is unavailable
    cg = (1/mPnts)*sum(pnts,2);

    [U,S,~] = svd( bsxfun(@minus,pnts,cg) );

    % Find minimum singular value
    s = diag(S);
    idx = find(s == min(s),1);

    n = U(:,idx);   % unit normal
    p = -n.'*cg;    % intercept
else
    % Fit nD plane to mPnts = nD if nCross is available
    cg = (1/mPnts)*sum(pnts,2);

    for i = 2:nD
        v{i-1} = pnts(:,i) - pnts(:,i-1);
    end
    n = nCross(v{:});

    n = n./norm(n);     % unit normal
    p = -n.'*pnts(:,1); % intercept
end
pln = zeros(1,numel(n)+1);
pln(1:numel(n)) = n;
pln(numel(n)+1) = p;

%% Package outputs
if nargout > 0
    varargout{1} = pln;
end

if nargout > 1
    [~,meanError] = proj2plane(pln,pnts);
    varargout{2} = meanError;
end