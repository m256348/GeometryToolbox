function varargout = fitCircle(X)
% FITCIRCLE fits a circle to a three-dimensional set of data.
%   cfit = FITCIRCLE(X) fits a circle to set of N three-dimensional points.
%   
%       X    - 3xN array containing points
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%
%   [..., meanError] = FITCIRCLE(...) additionally returns the mean error 
%   of the distance betweem provided points and their corresponding 
%   projections to the fit.
%
%   See also
%
%   M. Kutzer, 20Dec2017, USNA

%% Check Inputs
narginchk(1,1);

% Check points
[m,n] = size(X);
if m ~= 3
    error('Specified points must be provided as a 3xN array.');
end

%% Fit plane and define normal
pln = fitPlane(X);
n = transpose( pln(1:3) );
n_hat = n./norm(n);

%% Project points to plane
[Xproj,~,D] = proj2plane(pln,X);

%% Rotate points so they are parallel to the x/y plane
z_hat = n_hat;
x_hat = orthogonalVector(z_hat);
y = cross(z_hat,x_hat);
y_hat = y./norm(y);

R = [x_hat, y_hat, z_hat];

X0 = transpose(R)*Xproj;

%% Fit 2D circle
% Reference:
%   Izhak bucher, 25Oct1991
%   http://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit/content//circfit.m
x0 = transpose(X0(1,:));
y0 = transpose(X0(2,:));
a = [x0,y0,ones(size(x0))]\[-(x0.^2+y0.^2)];

% Define center in rotated frame
x0c = -0.5*a(1);
y0c = -0.5*a(2);
z0c = mean(X0(3,:));
C0 = [x0c; y0c; z0c];
% Define radius
r  =  sqrt((a(1)^2+a(2)^2)/4-a(3));

%% Define center relative to original coordinate system
C = R*C0;

%% Package outputs
cfit.Center = C;
cfit.Normal = n_hat;
cfit.Radius = r;

if nargout > 0
    varargout{1} = cfit;
end

if nargout > 1   
    % Calculate error in plane
    v = Xproj - C;
    ri = sqrt( sum( v.^2,1) );
    di = ri - r;
    
    % Compine errors
    err = sqrt( D.^2 + di.^2 );
    meanError = mean(err);
    varargout{2} = meanError;
end
