function bin = inEllipsoid(efit,X,ZERO)
% INELLIPSOID returns a binary array indicating whether points are
% contained within the volume of the specified ellipsoid.
%   in = INELLIPSOID(efit,X) identifies all points contained within the
%   volume of a specified ellipsoid or on the surface with a binary value 
%   of 'true'.
%   
%       X    - 3xN array containing points
%       efit - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%       in  - 1xN binary array containing 'true' if the corresponding
%           point is contained within or on the ellipsoid
%
%
%   [in,on] = INELLIPSOID(___) also returns a 1xN binary array indicating 
%   what points are exclusively on the ellipsoid. 
%
%   [___] = INELLIPSOID(___,ZERO) identifies the acceptable floating point
%   threshold. If this parameter is not specified, eps.m is used to
%   estimate floating point relative accuracy.
%
%   See also fitEllipsoid, proj2ellipsoid, makeEllipsoid
%
%   M. Kutzer, 03Jan2018, USNA

%% Check inputs
narginchk(2,3);

% Check efit
flds = {'Center','Rotation','PrincipalRadii'};
if ~all( isfield(efit,flds) )
    msg = 'Ellipsoid must be specified by a structured array containing the following fields: ';
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

%% Translate/Rotate points to body-fixed frame of ellipsoid.
X0 = X - efit.Center;
X0 = transpose(efit.Rotation) * X0;

%% Find points within the ellipsoid
chk1 = ( X0(1,:)./efit.PrincipalRadii(1) ).^2 + ...
       ( X0(2,:)./efit.PrincipalRadii(2) ).^2 + ...
       ( X0(3,:)./efit.PrincipalRadii(3) ).^2;

if nargin < 3
    ZERO = max( 100*eps([chk1,1]) );
end

chk0 = chk1 - 1;
bin = chk0 < ZERO;