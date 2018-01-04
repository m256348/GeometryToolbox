function efit_b = transformEllipsoid(efit_a,H_a2b)
% TRANSFORMELLIPSOID transforms a designated ellipsoid to a new frame of
% reference.
%   efit_b = transformEllipsoid(efit_a,H_a2b) transforms an ellipsoid
%   referenced to frame a to the same ellipsoid referenced to frame b.
%
%       H_a2b - 4x4 rigid body transformation defining frame a relative to
%           frame b.
%       efit  - structured array containing the following fields
%           efit.Center         - 3x1 center of the ellipsoid
%           efit.Rotation       - 3x3 rotation of the ellipsoid
%           efit.PrincipalRadii - radii of each principal semi-axis
%
%   See also fitEllipsoid, proj2ellipsoid, makeEllipsoid, inEllipsoid
%
%   M. Kutzer, 03Jan2018, USNA

%% Check inputs
narginchk(2,2);

% Check efit
flds = {'Center','Rotation','PrincipalRadii'};
if ~(sum( isfield(efit_a,flds) ) == numel(flds))
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

% Check H_a2b
if ~isSE(H_a2b)
    error('Transformation must be a valid element of the special Euclidean group.');
end

%% Transform ellipsoid
R_o2a = efit_a.Rotation;
T_o2a = efit_a.Center;

H_o2a = eye(4);
H_o2a(1:3,1:3) = R_o2a;
H_o2a(1:3,4) = T_o2a;

H_o2b = H_a2b*H_o2a;

efit_b.Rotation = H_o2b(1:3,1:3);
efit_b.Center = H_o2b(1:3,4);
efit_b.PrincipalRadii = efit_a.PrincipalRadii;