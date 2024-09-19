function Xint = intersectLineSegment(abc,seg,ZERO)
% INTERSECTLINESEGMENT calculates the intersection between a line and
% segment
%   Xint = intersectLineSegment(abc,seg)
%
%   Input(s)
%       abc - 1x(N+1) array defining the coefficients to ND line
%       seg - Nx2 array containing coefficients for the ND segment
%      ZERO - [OPTIONAL] scalar value that is sufficiently close to zero.
%             The default value is 1e-8. 
%
%   Output(s)
%       Xint - Nx1 array defining the point of intersection if it exists.
%              * If the segment and line do not intersect, Xint = []. 
%              * If the segment and line overlap, Xint is a Nx2 defining
%                the end-points of the segment. 
%
%   See also fitSegment fitLine
%
%   M. Kutzer, 19Sep2024, USNA

%% Check input(s)
narginchk(2,3);

nLine = numel(abc);
nSeg = size(seg,1);
mSeg = size(seg,2);

if mSeg ~= 2
    error('Segment coefficients must be defined as an Mx2 array.');
end

if nLine-1 ~= nSeg
    error('For a line with coefficients defined using a %d-element array, the segment coefficients must be defined as a %dx2 array.',...
        nLine,nLine-1);
end

if nargin < 3
    ZERO = 1e-8;
else
    if numel(ZERO) ~= 1
        error('ZERO must be defined as a positive scalar.');
    end
    if ZERO < 0
        error('ZERO must be defined as a positive scalar.');
    end
end

%% Parse input(s)
abc = reshape(abc,1,[]);
ab = abc(1:end-1);
c = abc(end);

%% Check for parallel lines
if abs(ab*seg(:,1)) < ZERO
    warning('Line and segment are near parallel.');
    
    % Define 2-point candidate intersection
    Xint = seg*[0,1;1,1];
    
    % Check if intersection is on the line
    if abs( abc*[Xint;1,1] ) > ZERO
        Xint = [];
    else
        warning('Line and segment are near parallel and overlapping. Returning two points of intersection.');
    end
    
    return
end

%% Solve
% -(c + a*c12 + b*c22)/(a*c11 + b*c21)
s = -(c + ab*seg(:,2))/(ab*seg(:,1));

if s >= -ZERO && s <= 1+ZERO
    Xint = seg*[s;1];
else
    Xint = [];
end