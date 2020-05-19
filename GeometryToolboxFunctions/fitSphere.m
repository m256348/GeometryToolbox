function varargout = fitSphere(X)
% FITSPHERE fits a sphere to a three-dimensional set of data.
%   sfit = FITSPHERE(X) fits a sphere to a set of N three-dimensional
%   points.
%
%       X    - 3xN array containing points
%       sfit - a sphereModel object if the sphereModel class is available 
%           *or* a 1x4 array containing sphere parameters [a,b,c,r] such  
%           that (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2
%
%   [..., meanError] = FITSPHERE(...) additionally returns the mean error 
%   of the distance of provided points and fit.
%
%   References
%       [1] https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
%
%   See also proj2sphere, pcfitsphere, sphereModel
%
%   M. Kutzer, 21Dec2017, USNA

%% Check Inputs
narginchk(1,1);

% Check points
[m,n] = size(X);
if m ~= 3
    error('Specified points must be provided as a 3xN array.');
end

%% Center data at the mean
X_bar = mean(X,2);
X0 = X - X_bar;

%% Estimate coefficients to quadratic equation
% See 7.2 in [1]
r2 = sum(X0.^2, 1);
x = X0(1,:);
y = X0(2,:);
z = X0(3,:);

M_hat(1,:) = [      n,         0,          0,          0,  sum(r2)];
M_hat(2,:) = [      0, sum(x.*x),  sum(x.*y),  sum(x.*z),  sum(x.*r2)];
M_hat(3,:) = [      0, sum(y.*x),  sum(y.*y),  sum(y.*z),  sum(y.*r2)];
M_hat(4,:) = [      0, sum(z.*x),  sum(z.*y),  sum(z.*z),  sum(z.*r2)];
M_hat(5,:) = [sum(r2), sum(r2.*x), sum(r2.*y), sum(r2.*z), sum(r2.*r2)];

%% Compute C
% See 7.2 in [1]
% C corresponds to the unit-length eigenvector corresponding to the minimum
% eigenvalue
[V,~] = eig(M_hat);
C = V(:,1);

%% Computer Center
% See 7.2 in [1]
% NOTE: C = [c0,c1,c2,c3,c4]
Center = transpose( X_bar - C(2:4)./(2*C(5)) );
Radius = abs( sqrt( sum(C(2:4).^2) - 4*C(1)*C(5))/(2*C(5)) );

%% Package outputs
sfit = [Center, Radius];

if nargout > 0
    if exist('sphereModel','class') == 8
        varargout{1} = sphereModel(sfit);
    else
        varargout{1} = sfit;
    end
end

if nargout > 1   
    % Calculate mean error
    [~,meanError] = proj2sphere(sfit,X);
    varargout{2} = meanError;
end