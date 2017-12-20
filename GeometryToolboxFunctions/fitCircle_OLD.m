function varargout = fitCircle(X)
% FITCIRCLE fits a circle to a three-dimensional set of data.
%   cfit = fitCircle(X) fits a circle to set of N three-dimensional points.
%   
%       X    - 3xN array containing points
%       cfit - structured array containint the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%
%   [..., meanError] = additionally returns the mean error of the distance
%   of provided points and fit.
%
%   References
%       [1] D. Elberly, "Fitting 3D Data with a Cylinder," 
%           https://www.geometrictools.com/Documentation/CylinderFitting.pdf
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

%% Subtract the sample mean
X_mean = mean(X,2);
X0 = X - X_mean;

%% Define normal
pln = fitPlane(X0);
W = transpose( pln(1:3) );
W = W./norm(W);

%% Solve for center
I = eye(3);
P = I - W*transpose(W);
S = wedge(W);

A = zeros(3,3);
B = zeros(3,1);
for i = 1:n
    Xi = X(:,i);
    A = A + Xi*transpose(Xi);
    B = B + (transpose(Xi)*P*Xi)*P*Xi;
end
A = P*A;

A_hat = S*A*transpose(S);

PC = (A_hat/trace(A_hat*A))*B;

PC
C = PC + X_mean;

%% Solve for radius
ri2 = zeros(1,n);
for i = 1:n
    Xi = X(:,i);
    Xi
    C
    P
    ri2(i) = transpose(C-Xi)*P*(C-Xi);
end

r2 = mean(ri2);
r = sqrt(r2);

%% Package outputs
cfit.Center = C;
cfit.Normal = W;
cfit.Radius = r;

if nargout > 0
    varargout{1} = cfit;
end

if nargout > 1
    % Project points and calculate error to plane of circle
    d = -dot(W,C);
    pln = [transpose(W),d];
    [Xproj,~,D] = proj2plane(pln,X);
    
    % Calculate error in plane
    v = Xproj - C;
    ri = sqrt( sum( v.^2,1) );
    di = ri - r;
    
    % Compine errors
    err = sqrt( D.^2 + di.^2 );
    meanError = mean(err);
    varargout{2} = meanError;
end