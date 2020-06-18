function varargout = inSphere(sfit,X)
% INSPHERE returns a binary array indicating whether points are
% contained within the volume of the specified sphere.
%   in = INSPHERE(sfit,X) returns a 1xN binary array indicating what
%   points are in or on the sphere. 
%
%       X    - 3xN array containing points
%       sfit - a sphereModel object if the sphereModel class is available 
%           *or* a 1x4 array containing sphere parameters [a,b,c,r] such  
%           that (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2
%
%   [in,on] = INSPHERE(___) also returns a 1xN binary array indicating what
%   points are exclusively on the sphere. 
%
%   [___] = INSPHERE(___,ZERO) identifies the acceptable floating point
%   threshold. If this parameter is not specified, eps.m is used to
%   estimate floating point relative accuracy.
%
%   See also fitSphere, sphereParam2Model, proj2sphere
%
%   M. Kutzer, 18Jun2020, USNA

%% Check inputs
narginchk(2,2);

% Standardize input
sfit = sphereParam2Model(sfit);

if size(X,1) ~= 3
    error('Specified points must be provided as a 3xN array.');
end

%% Translate points to a body-fixed frame centered in the sphere
X0 = X - reshape(sfit.Center,[],1);

%% Find points within the sphere
chk1 = ( X0(1,:)./sfit.Radius ).^2 + ...
       ( X0(2,:)./sfit.Radius ).^2 + ...
       ( X0(3,:)./sfit.Radius ).^2;
   
if nargin < 3
    ZERO = max( 100*eps([chk1,1]) );
end

chk0 = chk1 - 1;
in = chk0 < ZERO;

if nargout > 0
    varargout{1} = in;
end
if nargout > 1
    on = abs(chk0) < ZERO;
    varargout{2} = on;
end
