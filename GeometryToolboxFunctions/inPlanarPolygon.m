function [in,on] = inPlanarPolygon(Xq,Xv)
% INPLANARPOLYGON determines whether 3D points are contained within a
% planar polygon defind in 3D. 
%   in = INPLANARPOLYGON(Xq,Xv) returns "in" indicating if the query points
%   specified by the Nx3 array Xq are inside or on the edge of the polygon
%   area defined by the Mx3 array Xv. 
%
%   [in,on] = INPLANARPOLYGON(Xq,Xv) also returns "on" indicating if the
%   query points are on the edge of the polygon area.
%
%   See also inpolygon
%
%   M. Kutzer, 12Jun2020, USNA

% TODO - define a smarter value for "zero"
ZERO = 1e-12;

%% Check inputs
narginchk(2,2);

%% Fit a plane to Xv
abcd = fitPlane(Xv.');

%% Project points to the plane
Xv = proj2plane(abcd,Xv.').';

%% Define polygon-fixed coordinate frame
% Body-fixed z-direction
z_hat = abcd(1:3);
% Body-fixed origin
x0 = Xv(2,:).';
% Body-fixed 
x_hat = Xv(1,:).' - x0;
x_hat = x_hat./norm(x_hat);
y_hat = cross(z_hat,x_hat);
% Combine transformation
H_p2w = eye(4);
H_p2w(1:3,:) = [x_hat; y_hat; z_hat; d].';

%% Transform points
H_w2p = invSE(H_p2w);

Xq_w = Xq.';
Xq_w(4,:) = 1;
Xq_p = H_w2p * Xq_w;

Xv_w = Xv.';
Xv_w(4,:) = 1;
Xv_p = H_w2p * Xv_w;

% Find points that are effectively in plane
inPlane = abs(Xq_p(3,:)) < ZERO;
% Find 2D points that are in/on the polygon
[in,on] = inpolygon(Xq_p(1,:),Xq_p(2,:),Xv_p(1,:),Xv_p(2,:));

%% Combine for output
in = (in & inPlane).';
on = (on & inPlane).';