function cfit = intersectPlaneSphere(abcd,sfit)
% INTERSECTPLANESPHERE finds the parameters of the circle that lies at the
% intersection of a plane and sphere (if it exists).
%   cfit = INTERSECTPLANESPHERE(abcd,sfit)
%
%   Inputs:
%       abcd - 1x4 array containing the coefficients of a plane
%       sfit - sphereModel object or structured array containing
%           sfit.Center - 1x3 center of the sphere
%           sfit.Radius - scalar radius of the sphere
%   Outputs:
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%
%       NOTE: cfit is an empty set if no intersection exists
%
%   M. Kutzer, 22Jun2020, USNA

%% Check inputs
narginchk(2,2);

if numel(abcd) ~= 4
    error('Plane must be defined using 4 coefficients.');
end
abcd = reshape(abcd,1,[]);

% Standardize input
sfit = sphereParam2Model(sfit);

%% Find center
cfit.Center = reshape(proj2plane(abcd,reshape(sfit.Center,[],1)),[],1);

%% Define normal
cfit.Normal = reshape(abcd(1:3)./norm(abcd(1:3)),[],1);

%% Find distance between centers
d = sqrt( sum( (sfit.Center - cfit.Center.').^2 ) );

if d > sfit.Radius
    % No intersection exists
    cfit = [];
    return;
end

%% Find radius
cfit.Radius = sqrt( sfit.Radius.^2 - d.^2 );