%% SOLVE_intersectLineSegment
% Solve for the intersection between a line and parameterized segment
%
%   M. Kutzer, 19Sep2024, USNA

syms c11 c12 c21 c22 a b c s

C = [c11 c12; c21 c22];
x = C(1,:)*[s;1];
y = C(2,:)*[s;1];

abc = [a,b,c];

solve(abc*[x;y;1] == 0,s)