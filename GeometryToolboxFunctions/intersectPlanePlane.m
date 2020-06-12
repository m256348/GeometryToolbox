function M = intersectPlanePlane(abcd1,abcd2)
% INTERSECTPLANEPLANE calculates the line representing the intersection
% between two planes.
%   M = INTERSECTPLANEPLANE(abcd1,abcd2) specifies two planes as 1x4 arrays
%   of coefficients such that a*x + b*y + c*z + d = 0. This function
%   returns a 3x2 array M defining the coefficients of a parametric line in
%   3D:
%       X = M*[s; 1];
%
%   References:
%       [1] https://mathworld.wolfram.com/Plane-PlaneIntersection.html
%       [2] Gellert, W.; Gottwald, S.; Hellwich, M.; Kästner, H.; and 
%           Künstner, H. (Eds.). VNR Concise Encyclopedia of Mathematics, 
%           2nd ed. New York: Van Nostrand Reinhold, pp. 541-543, 1989.
%
%   M. Kutzer, 12Jun2020, USNA

%% Check/parse inputs
narginchk(2,2);

if numel(abcd1) ~= 4 || numel(abcd2) ~= 4
    error('Planes must be specified using an array of four coefficients.');
end

abcd1 = reshape(abcd1,1,[]);
abcd2 = reshape(abcd2,1,[]);

%% Calculate intersection
% m*x0 = b
n1 = abcd1(1:3).';
n2 = abcd2(1:3).';
m = [n1,n2].';
b = -[abcd1(4); abcd2(4)];

v0 = null(m);
if size(v0,2) == 1
    M(:,1) = v0;        % Slope v0
    M(:,2) = pinv(m)*b; % Offset x0
else
    M = [];
end