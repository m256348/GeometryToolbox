function M = fitSegment(p1,p2)
% FITSEGMENT fits an n-dimensional parametric equation for a segment
% between two points.
%   M = FITSEGMENT(p1,p2) fits a segment between p1 and p2 such that:
%       p(s) = M*[s; 1] where s \in [0,1]
%       p(0) = M*[0; 1] = p1;
%       p(1) = M*[1; 1] = p2;
%
%   Inputs:
%       p1 - n-element array representing the start-point in n-dimensional 
%            space
%       p2 - n-element array representing the end-point in n-dimensional
%            space
%
%   Outputs:
%       M - nx2 array containing coefficients for the segment
%
%   M. Kutzer, 23Jul2020, USNA

%% Check inputs
if numel(p1) ~= numel(p2)
    error('Specified points must be the same dimension.');
end

% Make sure points are nx1
p1 = reshape(p1,[],1);
p2 = reshape(p2,[],1);

%% Fit segment
M(:,1) = p2 - p1;
M(:,2) = p1; 
