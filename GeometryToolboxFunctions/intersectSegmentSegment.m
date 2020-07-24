function [intEE,intEV,intVV,intPt] = intersectSegmentSegment(p1,p2)
%INTERSECTSEGMENTSEGMENT finds intersection(s) between two segments.
%   [intEE,intEV,intVV,intPt] = INTERSECTSEGMENTSEGMENT(p1,p2) fits a
%   parametric line segment from p1(:,1) to p1(:,2), and a segment from 
%   p2(:,1) to p2(:,2); calculates the intersection conditions; and 
%   calculates the point(s) of intersection when applicable.
%
%   Inputs:
%       p1 - nx2 array containing n-dimensional segment end-points as 
%            column vectors
%       p2 - nx2 array containing n-dimensional segment end-points as 
%            column vectors
%
%       In 2D:
%           p1 = [x1_1, x1_2; y1_1, y1_2]
%           p2 = [x2_1, x2_2; y2_1, y2_2]
%
%   Outputs:
%       intEE - scalar binary that is true if there is at least one 
%               edge/edge intersect
%       intEV - scalar binary that is true if there is at least one 
%               edge/vertex intersect
%       intVV - scalar binary that is true if there is at least one 
%               vertex/vertex intersect
%       intPt - array containing coordinate(s) of the point or points  
%               of intersection (empty set, nx1, or nx2 array).
%           intPt = [];                 % No intersection exists
%           intPt = [x; y; ...];        % One intersection exists *or* 
%                                       % extending the segments into lines
%                                       % yields one intersection (this is 
%                                       % for debugging)
%           intPt = [x1,x2;y1,y2; ...]; % Two intersections exist (the  
%                                       % segments are parallel, colinear,  
%                                       % and overlapping)
%
%   Output Combinations:
%       (1 ) intEE = 0; intEV = 0; intVV = 0;
%       (2 ) intEE = 1; intEV = 0; intVV = 0;
%       (3 ) intEE = 0; intEV = 1; intVV = 0;
%       (4 ) intEE = 0; intEV = 0; intVV = 1;
%       (5a) intEE = 1; intEV = 1; intVV = 0;
%       (5b) intEE = 1; intEV = 0; intVV = 1;
%            intEE = 0; intEV = 1; intVV = 1; <- NOT POSSIBLE
%       (5c) intEE = 1; intEV = 1; intVV = 1;
%
%   Output Combinations Explained:
%       (1) The segments do not intersect:
%           -> intEE = 0; intEV = 0; intVV = 0;
%       (2) The edge of each segment intersects:
%           -> intEE = 1; intEV = 0; intVV = 0;
%       (3) The edge of one segment intersects with the vertex of the
%           other:
%           -> intEE = 0; intEV = 1; intVV = 0;
%       (4) A vertex of one segment is co-located with a vertex of the
%           other segment:
%           -> intEE = 0; intEV = 0; intVV = 1;
%       (5) The segments are parallel, co-linear, and overlapping:
%           (a) A vertex of one segment lies on the edge of the other 
%               segment and vice versa, or the both vertices of one segment 
%               lie on the other segment:
%               -> intEE = 1; intEV = 1; intVV = 0;
%           (b) The segments share both vertices:
%               -> intEE = 1; intEV = 0; intVV = 1;
%           (c) The segments share a vertex, and the vertex of one segment
%               lies on the edge of the other segment:
%               -> intEE = 1; intEV = 1; intVV = 1;
%
%   Example(s):
%       (1) Check if two segments intersect:
%       p1 = rand(2,2);
%       p2 = rand(2,2);
%       [intEE,intEV,intVV] = intersectSegmentSegment(p1,p2);
%       if any( [intEE,intEV,intVV] )
%           fprintf('Segments intersect.\n');
%       else
%           fprintf('Segments do not intersect.\n');
%       end
%
%   See also fitSegment
%
%   M. Kutzer, 18Mar2018, USNA

% Updates:
%   23Jul2020 - Full overhaul to address:
%               (1) singularity issues and 
%               (2) special case conditions associated with parallel lines
%   24Jul2020 - Extended to n-dimensional segments

%% Define "zero"
ZERO = 1e-6;    % TODO - use zeroFPError.m

%% Initialize defaults
intEE = false;  % No edge/edge intersections
intEV = false;  % No edge/vertex intersections
intVV = false;  % No vertex/vertex intresections
intPt = [];     % No points of intersection

%% Find the coefficients for the segments
M1 = fitSegment(p1(:,1),p1(:,2));
M2 = fitSegment(p2(:,1),p2(:,2));

%% Isolate slope/offset terms
A1 = M1(:,1); % Slope of segment 1
B1 = M1(:,2); % Offset of segment 1
A2 = M2(:,1); % Slope of segment 2
B2 = M2(:,2); % Offset of segment 2

%% Define "slope" and "offset" matrices for finding the intersection
%   A1*s1 + B1 = A2*s2 + B2
%   A1*s1 - A2*s2 = B2 - B1
%   [A1,-A2]*[s1; s2] = (B2 - B1)
%   AA*[s1; s2] = BB
AA = [A1, -A2]; % "slope" matrix
BB = (B2 - B1); % "offset matrix

%% Special Case: Segments are Parallel
% Check if lines are parallel

%if abs( det(AA) ) < ZERO % <- Only works for 2D points

% Segments are parallel if:
%   dot(A1,A2) = +/- norm(A1)*norm(A2)
if abs( abs(dot(A1,A2)) - abs(norm(A1)*norm(A2)) ) < ZERO
    % Lines are parallel!
    %warning('Segments are parallel or near parallel. det(AA) = %.10f',det(AA));
    
    % Check if segments overlap
    M3 = fitSegment(B1,B2); % Fit a segment to the "offsets" of each segment
    A3 = M3(:,1);           % Slope of "offsets" segment
    % Segments are colinear if (but we only need to check one):
    %   dot(A1,A3) = +/- norm(A1)*norm(A3)
    %   dot(A2,A3) = +/- norm(A2)*norm(A3)
    if abs( abs(dot(A1,A3)) - abs(norm(A1)*norm(A3)) ) < ZERO
        % Segments are colinear!
        
        % One of the following will apply:
        %   (1) Segments do not overlap
        %   (2) Segments overlap at exactly one point
        %   (3) Segments overlap over a range of points
        
        % Define points of intersection
        nZ = (A1 ~= 0); % Consider non-zero values of A1 (avoid divide by 0)
        S1 = (p2(nZ,:) - repmat(B1(nZ,:),1,2))./A1(nZ,:); % End-points of segment 2 lie on segment 1
        nZ = (A2 ~= 0); % Consider non-zero values of A2 (avoid divide by 0)
        S2 = (p1(nZ,:) - repmat(B2(nZ,:),1,2))./A2(nZ,:); % End-points of segment 1 lie on segment 2
        % Each element of a column *should* be identical
        S1 = mean(S1,1);
        S2 = mean(S2,1);
        
        p2On1 = ( (S1+ZERO) >= 0 & (S1-ZERO) <= 1); % Segment 2 end-points that lie on Segment 1
        p1On2 = ( (S2+ZERO) >= 0 & (S2-ZERO) <= 1); % Segment 1 end-points that lie on Segment 2
        
        % (1) Segments do not overlap
        if all(~p2On1) && all(~p1On2)
            % (1) Segments do not overlap
            intEE = false;  % No edge/edge intersects
            intEV = false;  % No edge/vertex intersects
            intVV = false;  % No vertex/vertex intersects
            intPt = [];     % Return no points
            return
        end
        
        % (2) Segments overlap at exactly one point
        % (3) Segments overlap over a range of points
        intPt = [p1(:,p1On2),p2(:,p2On1)]; % There should be between 2 and 4 points returned
        switch size(intPt,2)
            case 4
                % Segments share both vertices
                % (3) Segments overlap over a range of points
                %     -> Entire segment is shared
                intEE = true;  % Infinite edge/edge intersects
                intEV = false; % No edge/vertex intersects
                intVV = true;  % Two vertex/vertex intersect
                intPt = p1;    % Return two points
                return
            case 3
                % Segments share one vertex and one vertex lies on an edge
                % (3) Segments overlap over a range of points
                intEE = true;   % Infinite edge/edge intersects
                intEV = true;   % One edge/vertex intersect
                intVV = true;   % One vertex/vertex intersect
                if nnz(p1On2) == 2
                    intPt = p1; % Return two points
                else
                    intPt = p2; % Return two points
                end
                return
            case 2
                % Two vertices lie on an edge OR
                % Segments share one single vertex
                if norm(diff(intPt,1,2)) < ZERO
                    % (2) Segments overlap at exactly one point
                    intEE = false;      % No edge/edge intersects
                    intEV = false;      % No edge/vertex intersects
                    intVV = true;       % One vertex/vertex intersect
                    intPt = intPt(:,1); % Return one point
                    return
                else
                    % (3) Segments overlap over a range of points
                    intEE = true;  % Infinite edge/edge intersects
                    intEV = true;  % Two edge/vertex intersects
                    intVV = false; % No vertex/vertex intersects
                    intPt = intPt; % Return two points
                    return
                end
            otherwise
                % Dump variables
                p1, M2, S2, p1On2, ...
                p2, M1, S1, p2On1, ... 
                intPt
                % Throw error
                error('Unexpected number of intersect points found.');
        end
    else
        % Segments are parallel but not colinear
        %   -> No intersection exists
        intEE = false;  % No edge/edge intersections
        intEV = false;  % No edge/vertex intersections
        intVV = false;  % No vertex/vertex intresections
        intPt = [];     % No points of intersection
        return
    end
end

%% Special Case: Non-parallel vertex/vertex intersection
% Enumerate all combinations of p1 and p2
enum_p1 = [p1(:,1),p1(:,2),p1(:,1),p1(:,2)];
enum_p2 = [p2(:,1),p2(:,1),p2(:,2),p2(:,2)];
% Fund the sum square difference of all combinations
enum_ssd = sum((enum_p1 - enum_p2).^2, 1); 
if any( enum_ssd < ZERO )
    intEE = false; % No edge/edge intersections
    intEV = false; % No edge/vertex intersections
    intVV = true;  % one vertex/vertex intresections
    intPt = enum_p1(:,enum_ssd < ZERO); % One point of intersection
    return
end

%% Calculate s-values of intersection
% Find parametric intersection
if size(AA,1) == size(AA,2)
    % Two-dimensional Segment
    s1s2 = AA^(-1) * BB;
else
    % n-dimensional segment (n > 2)
    s1s2 = pinv(AA) * BB;
    % Check "intersection" to see if it is the same point
    if norm( M1*[s1s2(1); 1] - M2*[s1s2(2); 1] ) > ZERO
        % No actual intersection exists
        intEE = false; % No edge/edge intersections
        intEV = false; % No edge/vertex intersections
        intVV = false; % No vertex/vertex intresections
        intPt = [];    % No points of intersection
        return
    end
end

%% Check parameter bounds 
IsOn = ( (s1s2+ZERO) >= 0 & (s1s2-ZERO) <= 1); % Intersections are on the segment
NotOn = ~IsOn;                                 % Intersections are *not* on the segment

%% Check if intersection occurs off of at least one segment (i.e. no actual intersection)
if any(NotOn)
    %fprintf('Intersect off-segment\n')
    % Intersection is not on at least one of the segments
    intEE = false; % No edge/edge intersections
    intEV = false; % No edge/vertex intersections
    intVV = false; % No vertex/vertex intresections
    intPt = M1*[s1s2(1);1]; % No point of intersection (return the line intersection for debugging)
    return
end

%% Check for edge/vertex intersection
OnEdge = ( abs(s1s2) < ZERO | abs(s1s2-1) < ZERO );
if any(OnEdge)
    if nnz(OnEdge) == 2
        % Dump variables
        p1, M2,...
        p2, M1,...
        s1s2
        % Throw warning
        warning('Vertex/vertex intersect found after special case!');
        intEE = false; % No edge/edge intersections
        intEV = false; % No edge/vertex intersections
        intVV = true;  % One vertex/vertex intresections
        intPt = M1*[s1s2(1);1]; % One point of intersection
        return
    else
        % Edge/vertex intersection
        intEE = false; % No edge/edge intersections
        intEV = true;  % One edge/vertex intersections
        intVV = false; % No vertex/vertex intresections
        intPt = M1*[s1s2(1);1]; % One point of intersection
        return
    end
end

%% Edge/edge intersection
% If we have made it this far into the code, the intersection must be
% edge/edge

% Edge/edge intersection
intEE = true;       % One edge/edge intersection
intEV = false;      % No edge/vertex intersections
intVV = false;      % No vertex/vertex intresections
intPt = M1*[s1s2(1);1]; % One point of intersection