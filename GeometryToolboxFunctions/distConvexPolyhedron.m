function [d,msg] = distConvexPolyhedron(A,B,ZERO)
% DISTCONVEXPOLYHEDRON calculates the distance between two convex
% n-dimensional polyhedron.
%
% !!!!!!!!!! THIS FUNCTION IS INCOMPLETE !!!!!!!!!!!
%
%   [d,msg] = DISTCONVEXPOLYHEDRON(A,B)
%   [d,msg] = DISTCONVEXPOLYHEDRON(A,B,ZERO)
%
%   Input(s)
%       A    - structured array (or object) containing fields "Faces" and 
%              "Vertices" or a valid array of patch objects.
%       B    - structured array containing fields "Faces" and "Vertices" 
%              or a valid array of patch objects.
%       ZERO - [OPTIONAL] value sufficiently close to zero to consider 
%              values the same. If unspecified, 0 is used.
%
%       *.Vertices - NxM array containing N, M-dimensional vertices
%       *.Faces    - nxm array containing n faces (defined using vertex
%                    index values)
%
%           Faces can consist of varying numbers of vertices.
%           E.g. Face 2 consists of two vertices, Face 3 consists of 5
%                vertices:
%               *.Faces(3,:) = [___, ___, ___, nan, nan]
%               *.Faces(5,:) = [___, ___, ___, ___, ___]
%
%   Output(s)
%       d   - scalar value describing distance. d = 0 if the polyhedron
%           overlap, are adjacent, or a single-point contact exists.
%       msg - string describing distance result:
%
%       Values of msg:
%       (1) Condition: d > 0
%           'Non-overlapping' - polyhedron are non-overlapping
%       (2) Condition: d <= 0 (0 is returned)
%               'SinglePoint' - polyhedron share a single vertex
%           'SinglePointAonB' - single vertex of polyhedron A lies on B
%           'SinglePointBonA' - single vertex of polyhedron B lies on A
%             'AdjacentFaces' - single face adjacency of polyhedron A/B 
%             'AdjacentEdges' - single edge adjacency of polyhedron A/B 
%                   'Overlap' - polyhedron A and B overlap
%                      'AinB' - polyhedron A lies entirely inside B
%                      'BinA' - polyhedron B lies entirely inside A
%                    'AinonB' - polyhedron A lies inside and on B
%                    'BinonA' - polyhedron B lies inside and on A
%
%   M. Kutzer, 06Apr2022, USNA

%% Check input(s)
% NOTE: We are not going to check the following:
%   (1) Faces containing more vertices than necessary to define the
%       hyperplane of the face, we will not check to confirm that all
%       points are co-planar
%   (2) Face ordering to ensure that vertices produce an outward pointing
%       normal
%   (3) Whether the polyhedron is convex
%   (4) Whether the polyhedron is closed
%
% If any of these properties are violated, the results will be incorrect.
narginchk(2,3);

% Set default for ZERO
if nargin < 3
    ZERO = [];
end

% Check for and isolate faces and vertices
try
    vertsA = A.Vertices;
    vertsB = B.Vertices;
    facesA = A.Faces;
    facesB = B.Faces;
catch
    error('A and B must contain the fields or properties "Faces" and "Vertices".');
end

% Check dimensions 
nV_A = size(vertsA,1);  % Number of vertices, Polyhedron A
nD_A = size(vertsA,2);  % Dimension of vertices, Polyhedron A
nV_B = size(vertsB,1);  % Number of vertices, Polyhedron B
nD_B = size(vertsB,2);  % Dimension of vertices, Polyhedron A
if nD_A ~= nD_B
    error(...
        ['Vertices of A and B must be the same dimension\n',...
        '\tVertices of A are %d dimensional\n',...
        '\tVertices of B are %d dimensional'],nD_A,nD_B);
end
nD = nD_A;  % Common dimension of vertices

% Check for minimum number of faces to define hyperplanes
% TODO - check for minimum number of faces

%% Apply ZERO
if ~isempty(ZERO)
    vertsA = round(vertsA./ZERO).*ZERO;
    vertsB = round(vertsB./ZERO).*ZERO;
end

%% Fit hyperplanes to each face
nFaces_A = size(facesA,1);
nFaces_B = size(facesB,1);

facePlaneA = nan(nFaces_A,nD+1);
facePlaneB = nan(nFaces_B,nD+1);
for i = 1:max([nFaces_A,nFaces_B])
    % Face planes, Polyhedron A
    if i <= nFaces_A
        for iA = 1:nD
            pA(:,iA) = vertsA(facesA(i,iA),:).';
        end
        facePlaneA(i,:) = fitPlane(pA);
    end
    
    % Face planes, Polyhedron B
    if i <= nFaces_B
        for iB = 1:nD
            pB(:,iB) = vertsB(facesB(i,iB),:).';
        end
        facePlaneB(i,:) = fitPlane(pB);
    end
end

%% Check for vertices on the interior
xA = vertsA.';
xA(4,:) = 1;
xB = vertsB.';
xB(4,:) = 1;

% Check for values on the interior
AinB = (facePlaneB*xA) < 0;
BinA = (facePlaneA*xB) < 0;

% Check for values on the polygon
AonB = (facePlaneB*xA) == 0;
BonA = (facePlaneA*xB) == 0;

%% Categorize
if nnz(AinB) == nV_A
    msg = 'AinB';
    d = 0;
    return
end

if nnz(BinA) == nV_B
    msg = 'BinA';
    d = 0;
    return
end

if nnz(AinB)+nnz(AonB) >= nV_A
    msg = 'AinonB';
    d = 0;
    return
end

if nnz(BinA)+nnz(BonA) >= nV_B
    msg = 'BinonA';
    d = 0;
    return
end

if nnz(AinB) > 0 || nnz(BinA) > 0
    msg = 'Overlap';
    d = 0;
    return
end

% SinglePoint
% -> This assumes a single point adjacency is nD faces coinciding
if nnz(AonB) == nD && nnz(BonA) == nD
    msg = 'SinglePoint';
    d = 0;
    return
end

% TODO - Finish!
%   'SinglePointAonB'
%   'SinglePointBonA'
%   'AdjacentFaces'
%   'Non-overlapping'
%warning('THIS FUNCTION IS INCOMPLETE');
msg = 'INCOMPLETE!';
d = inf;