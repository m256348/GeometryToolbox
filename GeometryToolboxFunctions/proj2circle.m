function varargout = proj2circle(cfit,X)
% PROJ2CIRCLE projects a set of points to the edge of a specified circle
%   Xproj = PROJ2CIRCLE(cfit,X) projects the N points contained in the 3xN
%   array X to the edge of a circle specified as a structured array.
%   
%       X    - 3xN array containing points
%       cfit - structured array containing the following fields
%           cfit.Center - 3x1 center of the circle
%           cfit.Normal - 3x1 normal to the circle
%           cfit.Radius - radius of the circle
%
%   [..., meanError] = PROJ2CIRCLE(...) additionally returns the mean error
%   of the Euclidean distance between the provided points and their
%   corresponding projections.
%
%   [..., Error] = PROJ2CIRCLE(...) additionally returns the errors for
%   each projected point.
%
%   See also fitCircle
%
%   M. Kutzer, 03Jan2018, USNA

warning('THIS FUNCTION IS INCOMPLETE');

%% Check Inputs
narginchk(2,2);

% Check cfit
flds = {'Center','Normal','Radius'};
if ~(sum( isfield(cfit,flds) ) == numel(flds))
    msg = 'Circle must be specified by a structured array containing the following fields: ';
    for i = 1:numel(flds)
        if i < numel(flds)
            msg = sprintf('%s%s, ',msg,flds{i});
        else
            msg = sprintf('%sand %s.',msg,flds{i});
        end
    end
    error(msg);
end

% Check points
[M,~] = size(X);
if M ~= 3
    error('Specified points must be provided as a 3xN array.');
end