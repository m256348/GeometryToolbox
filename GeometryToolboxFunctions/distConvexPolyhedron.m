function [d,msg] = distConvexPolyhedron(A,B)
% DISTCONVEXPOLYHEDRON calculates the distance between two convex
% n-dimensional polyhedron.
%   [d,msg] = DISTCONVEXPOLYHEDRON(A,B)
%
%   Input(s)
%       A - structured array containing fields "Faces" and "Vertices" or a 
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
    
warning('THIS FUNCTION IS INCOMPLETE');