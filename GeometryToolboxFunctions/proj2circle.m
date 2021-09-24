function varargout = proj2circle(cfit,X)
% PROJ2CIRCLE projects a set of points to the edge of a specified circle
%   Xproj = PROJ2CIRCLE(cfit,X) projects the N points contained in the 3xN
%   array X to the edge of a circle specified as a structured array.
%   
%       X    - 3xN array containing points
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%
%   [..., meanError] = PROJ2CIRCLE(...) additionally returns the mean error
%   of the Euclidean distance between the provided points and their
%   corresponding projections.
%
%   [..., Error] = PROJ2CIRCLE(...) additionally returns the errors for
%   each projected point.
%
%   See also fitCircle
%
%   M. Kutzer, 03Jan2018, USNA

% Updates
%   23Sep2021 - Completed function

%% Check Inputs
narginchk(2,2);

% Check cfit
flds = {'Center','Normal','Radius'};
if ~(sum( isfield(cfit,flds) ) == numel(flds))
    msg = 'Circle must be specified by a structured array containing the following fields: ';
    for i = 1:numel(flds)
        if i < numel(flds)
            msg = sprintf('%s%s, ',msg,flds{i});
        else
            msg = sprintf('%sand %s.',msg,flds{i});
        end
    end
    error(msg);
end

% Check points
[M,~] = size(X);
if M ~= 3
    error('Specified points must be provided as a 3xN array.');
end

%% Define rigid body transform
% z-direction
z_hat = reshape(cfit.Normal./norm(cfit.Normal),[],1);
% "other" direction
n_hat = z_hat([2,3,1]);
% x-direction
x_hat = cross(n_hat,z_hat);
x_hat = x_hat./norm(x_hat);
% y-direction
y_hat = cross(z_hat,x_hat);

H_c2w = eye(4);
H_c2w(1:3,1:3) = [x_hat,y_hat,z_hat];
H_c2w(1:3,4) = reshape(cfit.Center,[],1);

%% Project points to plane
abcd = [z_hat.',-z_hat.'*cfit.Center];
X_w = proj2plane(abcd,X);

%% Reference points to body-fixed frame
X_w(4,:) = 1; % Make homogeneous
X_c = invSE(H_c2w)*X_w;

%% Project points
thetas = atan2(X_c(2,:),X_c(1,:));

pX_c(1,:) = cfit.Radius*cos(thetas);
pX_c(2,:) = cfit.Radius*sin(thetas);
pX_c(4,:) = 1;

%% Reference projected points to world frame
pX_w = H_c2w*pX_c;

%% Package outputs
varargout{1} = pX_w;
if nargout > 1
    dX = X_w(1:3,:) - pX_w(1:3,:);
    err = sqrt( sum( dX.^2,1) );
    meanError = mean(err);
    varargout{2} = meanError;
end
if nargout > 2
    varargout{3} = err;
end
    