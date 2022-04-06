function [d,msg] = distConvexPolyhedron(A,B)
% DISTCONVEXPOLYHEDRON calculates the distance between two convex
% n-dimensional polyhedron.
%   [d,msg] = DISTCONVEXPOLYHEDRON(A,B)
%
%   Input(s)
%       A - structured array (or object) containing fields "Faces" and "Vertices" or a 
%           valid array of patch objects.
%       B - structured array containing fields "Faces" and "Vertices" or a 
%           valid array of patch objects.
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
%           'Non-overlapping' - polygons are non-overlapping (d > 0)
%               'SinglePoint'
%           'SinglePointAonB'
%           'SinglePointBonA'
%             'AdjacentFaces'
%                   'Overlap'
%                      'AinB'
%                      'BinA'
%
%   M. Kutzer, 06Apr2022, USNA

warning('THIS FUNCTION IS INCOMPLETE');

%% Check input(s)
% NOTE: We are not going to check the following:
%   (1) For faces containing more vertices than necessary to define the
%       hyperplane of the face, we will not check to confirm that all
%       points are co-planar
%   (2) The face ordering ensuring that vertices produce an outward
%       pointing normal
%   (3) Whether the polyhedron is convex
%
% If any of these properties are violated, the results will be incorrect.
narginchk(2,2);

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
nD_A = size(vertsA,2);
nD_B = size(vertsB,2);
if nD_A ~= nD_B
    error(...
        ['Vertices of A and B must be the same dimension\n',...
        '\tVertices of A are %d dimensional\n',...
        '\tVertices of B are %d dimensional'],nD_A,nD_B);
end
nD = nD_A;

% Check for minimum number of faces to define hyperplanes
% TODO - check for minimum number of faces

%% Fit hyperplanes to each face
nFaces_A = size(facesA,1);
nFaces_B = size(facesB,1);

facePlaneA = nan(nFaces_A,nD+1);
facePlaneB = nan(nFaces_B,nD+1);
for i = 1:max([nFaces_A,nFaces_B])
    % Face planes, Polyhedron A
    if i <= nFacesA
        for iA = 1:nD
            pA(:,iA) = vertsA(facesA(i,iA),:).';
        end
        facePlaneA(i,:) = fitPlane(pA);
    end
    
    % Face planes, Polyhedron B
    if i <= nFacesB
        for iB = 1:nD
            pA(:,iB) = vertsB(facesB(i,iB),:).';
        end
        facePlaneB(i,:) = fitPlane(pB);
    end
end

%% Check half spaces
