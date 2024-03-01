function efit = mvnpdfConfInterval(mu,Sigma,prc)
% MVNPDFCONFINTERVAL defines the confidence interval for a multivariate 
% Gaussian distribution given mean, covariance, and a percent probability.
%
%   efit = mvnpdfConfInterval(mu,Sigma,prc)
%
%   Input(s)
%          mu - N-element array defining distribution mean
%       Sigma - NxN array defining distribution covariance
%         prc - scalar value defining probability associated with the  
%               confidence interval ( 0 > prc < 1 )
%
%   Output(s)
%       efit - structured array containing the following fields
%           efit.Center         - Nx1 center of the ellipsoid
%           efit.Rotation       - NxN rotation of the ellipsoid
%           efit.PrincipalRadii - Nx1 radii of each principal semi-axis
%
%   See also plotEllipsoid
%
%   M. Kutzer, 01Mar2024, USNA

%% Check input(s)
narginchk(3,3);

% Check distribution size
N = numel(mu);
if ~ismatrix(Sigma) || size(Sigma,1) ~= N || size(Sigma,2) ~= N
    error('Given a %d element mean, covariance must be a %dx%d matrix.',N,N,N);
end

% Check for positive definite covariance
try
    chol(Sigma);
catch
    error('Given covariance matrix must be symmetric positive definite.');
end

% Check percent confidence
if numel(prc) ~= 1 || prc <= 0 || prc >= 1
    error('Probability must be defined as a value between 0 and 1.');
end

%% Define value associated with confidence interval
val = sqrt( chi2inv(prc,N) );

%% Define Eigenvalues and Eigenvectors
[V,D] = eig(Sigma);

%% Define principal axis radii
r = val*sqrt( diag(D) );

%% Define rotation
if N == 1
    R_e2o = 1;
else
    % Define cross product combinations
    idxCross = flipud( nchoosek(1:N,N-1) );
    
    % Initialize unit vector directions
    for i = 1:N
        k{i} = V(:,i)./norm(V(:,i));
    end
    
    % Define Nth direction ensuring cross product
    for i = 2:N
        k{i} = nCross(k{idxCross(i,:)});
    end

    % Define rotation matrix
    R_e2o = eye(N);
    for i = 1:N
        R_e2o(:,i) = k{i};
    end
end

%% Package output(s)
efit.Center = reshape(mu,[],1);
efit.Rotation = R_e2o;
efit.PrincipalRadii = reshape(r,[],1);