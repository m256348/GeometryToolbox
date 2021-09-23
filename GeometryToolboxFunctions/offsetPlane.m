function pln1 = offsetPlane(pln,D)
% OFFSETPLANE defines a parallel plane offset by distance D.
%   pln1 = OFFSETPLANE(pln,D) defines a plane offset by distance D from a
%   specified plane.
%
%   Input(s)
%       pln - 1x(N+1) array containing coefficients for plane equation 
%              *or* a valid planeModel object
%          General: [c_1,c_2,c_3,... c_{N+1}] such that
%                      c_1*x + c_2*y + c_3*z ... + c_{N+1) = 0
%        Line (2D): [a,b,c] such that a*x + b*y + c = 0
%       Plane (3D): [a,b,c,d] such that a*x + b*y + c*z + d = 0 
%       D  - scalar offset
%
%   Output(s)
%       pln1 - 1x(N+1) array containing coefficients for plane equation *or* a valid
%              planeModel object 
%
%   M. Kutzer, 23Sep2021, USNA

%% Check Inputs
narginchk(2,2);

% Get plane normal
switch lower( class(pln) )
    case 'planemodel'
        n = transpose( pln.Normal );
        d = pln.Parameters(4);
    case 'double'
        n = transpose( pln(1:end-1) );
        d = pln(end);
    case 'single'
        n = transpose( pln(1:end-1) );
        d = pln(end);
    case 'sym'
        n = transpose( pln(1:end-1) );
        d = pln(end);
    otherwise
        error('Plane must be specified as a 1x(N+1) array where N>1 or a valid planeModel object.');
end

% Check normal
N = numel(n);
if N < 2
    error('Plane must be N-dimensional where N > 2.');
end

%% Define Hessian normal form of plane
%n_hat = n./norm(n);
p = d./norm(n);

%% Define updated p
p1 = p - D;

%% Package output
d1 = p1 * norm(n);

% Package plane
pln1 = [transpose(n),d1];
switch lower( class(pln) )
    case 'planemodel'
        pln1 = planeModel(pln1);
end