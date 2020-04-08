function [intEE,intEV,intVV,intPt] = intersectSegmentSegment(p1, p2)
%INTERSECTSEGMENTSEGMENT finds intersection(s) between two segments.
%   [intEE,intEV,intVV,intPt] = INTERSECTSEGMENTSEGMENT(p1, p2) calculates a
%   parametric line segment equation from p1(:,1) to p1(:,2) and p2(:,1) 
%   to p2(:,2); calculates the intersection conditions; and calculates the 
%   point(s) of intersection.
%
%   Inputs:
%       p1 = [x1_1, x1_2; y1_1, y1_2]
%       p2 = [x2_1, x2_2; y2_1, y2_2]
%
%   Outputs:
%       intEE - binary that is true if there is an edge/edge intersect
%       intEV - binary that is true if there is an edge/vertex intersect
%       intVV - binary that is true if there is an vertex/vertex intersect
%       intPt = [x; y] - the point of intersection (for debugging and stuff)
%
%   M. Kutzer, 18Mar2018, USNA

%% Initialize defaults
intEE = false;
intEV = false;
intVV = false;
intPt = [];

%% Fit our coefficients
M1 = p1*[0, 1; 1, 1]^(-1);
M2 = p2*[0, 1; 1, 1]^(-1);

%% Define "slope" matrix
MM = [M1(:,1), -M2(:,1)];

%% Define "offset" matrix
BB = (M2(:,2) - M1(:,2));

%% Check for parallel line condition
ZERO = 1e-6;    % TODO - use zeroFPError.m
isParallel = false;
if abs( det(MM) ) < ZERO
    % TODO - check if a row or column of M* is 0
    %   -> This corresponds to a vertical or horizontal segment that goes
    %   through the origin.
    
    % Lines are parallel!
    isParallel = true;
    
    % Check for overlapping lines
    for s = [0,1]
        xy1 = M1*[s; 1];  % Calculate end-point
        xy2 = M2*[s; 1];  % Calculate end-point
        
        S1 = M1^(-1) * xy2;
        S2 = M2^(-1) * xy1;
        
        if S1(1) > 0 && S1(1) < 1
            % Intersect
            intEE = true;
            intEV = true;
            intPt(:,end+1) = M1*S1;
        end
        if S2(1) > 0 && S2(1) < 1
            % Intersect
            intEE = true;
            intEV = true;
            intPt(:,end+1) = M2*S2;
        end
        % check for vertex/vertex (double check me)
        if (S1(1) == 1 || S1(1) == 0) && (S2(1) == 1 || S2(1) == 0)
            intVV = true;
            intPt(:,end+1) = M1*S1;
        end
    end
    
    % TODO - only return the unique set for intPt
    return
end
        
%% Calculate s-values
s1s2 = MM^(-1) * BB;
s1 = s1s2(1);
s2 = s1s2(2);

%% Check conditions
if s1 > 0 && s1 < 1 && s2 >0 && s2 < 1
    % Edge/Edge
    intEE = true;
end

if (s1 == 0 || s1 == 1) && (s2 == 0 || s2 == 1)
    % Vertex/Vertex
    intVV = true;
end

if (s1 == 0 || s1 == 1) && (s2 >0 && s2 < 1)
    % Vertex/Edge
    intEV = true;
end

if (s1 > 0 && s1 < 1) && (s2 == 0 || s2 == 1)
    % Edge/Vertex
    intEV = true;
end

%% Calculate the point of intersect
intPt = M1*[s1; 1];
%intPt = M2*[s2; 1];