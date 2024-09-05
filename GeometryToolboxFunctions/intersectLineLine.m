function X = intersectLineLine(abc1,abc2)
% INTERSECTLINELINE calculates the intersection between two lines
%   X = intersectLineLine(abc1,abc2)
%
%   Input(s)
%    abc1 - 1x3 array defining the coefficients to line 1
%    abc2 - 1x3 array defining the coefficients to line 2
%
%   Output(s)
%       X - 2x1 array defining point of intersection
%
%   See also fitPlane plotLine
%
%   M. Kutzer, 14May2024, USNA

ZERO = 1e-8;
%% Check input(s)
narginchk(2,2);

if numel(abc1) ~= 3 || ~isnumeric(abc1)
    error('Line 1 must be defined using three constants.');
end

if numel(abc2) ~= 3 || ~isnumeric(abc1)
    error('Line 1 must be defined using three constants.');
end

abc1 = reshape(abc1,1,[]);
abc2 = reshape(abc2,1,[]);

%% Find point of intersection
AB = [abc1(1:2); abc2(1:2)];
C = [abc1(3); abc2(3)];

if abs( det(AB) ) < ZERO
    warning('Lines are near parallel.');
end

X = -(AB^-1)*C;