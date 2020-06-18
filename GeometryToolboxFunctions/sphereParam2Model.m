function sfit = sphereParam2Model(params)
% SPHEREPARAM2MODEL converts various representations of sphere paramters to
% a sphere model or its structured array equivalent.
%   sfit = sphereParam2Model(params) returns a sphereModel object (if
%   available) or a structured array.
%
%       params - a sphereModel object various representations of a sphere 
%       including:
%           (1) a 1x1 array containing the radius of the sphere only; 
%           (2) a 1x4 array containing sphere parameters [a,b,c,r] such
%           that (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2; OR
%           (3) a structured array containing the fields "Center" and
%           "Radius".
%
%       sfit - a sphereModel object if the sphereModel class is available 
%           *or* a structured array containing the fields "Parameters," 
%           "Center," and "Radius."
%
%   M. Kutzer, 18Jun2020, USNA

%% Check inputs
narginchk(1,1);

switch class(params)
    case 'sphereModel'
        sfit = params;
        return
    case 'struct'
        strct = params;
        if ~all( isfield(strct,{'Radius','Center'}) )
            error('Structured array inputs must contain the fields "Radius" and "Center"');
        end
        params = [reshape(strct.Center,1,[]),strct.Radius];
end

%% Check for radius only
if numel(params) == 1
    % Assume [0,0,0] center
    params(1,4) = params;
end

%% Check params
if ~numel(params) == 4
    error('Parameters array must be 1x4');
end

params = reshape(params,1,[]);

if exist('sphereModel','class') == 8
    sfit = sphereModel(params);
else
    sfit.Parameters = params;
    sfit.Center = params(1:3);
    sfit.Radius = params(4);
end

