function varargout = plotSphere(varargin)
% PLOTSPHERE creates a patch object representation of a sphere
%   [...] = PLOTSPHERE(sfit)
%
%   [...] = PLOTSPHERE(axs,sfit) specifies the parent as axs
%
%   [...] = PLOTSPHERE(...,N) specifies the number of faces as N^2
%
%   p_a = PLOTSPHERE(...) creates a patch object representation of the 
%   sphere referenced to the default or specified parent (gca or axs)
%
%   [p_s,h_s2a] = PLOTSPHERE(...) creates a patch object representation of
%   the sphere centered at the origin parent transform h_s2a referenced to
%   the default or specified parent (gca or axs)
%
%   Input(s)
%       axs  - [OPTIONAL] parent object for the patch object representation
%              of the sphere. This can be an axes or hgtransform object. If
%              none is specified, "axs = gca" is used
%       sfit - a sphereModel object *or* 1x4 array containing sphere 
%              parameters [a,b,c,r] such that 
%              (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2
%
%           shereModel Object
%               sfit.Parameters = [a, b, c, r]
%               sfit.Center     = [a, b, c]
%               sfit.Radius     = r
%
%           Sphere Parameters
%               sfit = [a, b, c, r]
%
%       N    - [OPTIONAL] creates a sphere with N^2 faces (default N = 50)
%
%   Output(s)
%       p_a   - patch object referenced to the default or specified parent
%       p_s   - patch object centered at [0; 0; 0] and referenced to h_s2a
%       h_s2a - hgtransform relating p_s to the default or specified parent
%
%   M. Kutzer, 07Mar2022, USNA

% Updates
%   17Mar2022 - Corrected parent assignment error for non-axes parents

%% Check inputs
narginchk(1,3);

% Parse inputs
N = 50;
if nargin == 1
    axs = gca;
    sfit = varargin{1};
end
if nargin == 2
    if ishandle(varargin{1})
        axs = varargin{1};
        sfit = varargin{2};
    else
        axs = gca;
        sfit = varargin{1};
        N = varargin{2};
    end
end
if nargin == 3
    axs = varargin{1};
    sfit = varargin{2};
    N = varargin{3};
end

% Get plane normal
msg = [];
msgErr = 'Sphere must be specified as a 1x4 array or a valid sphereModel object.';
switch lower( class(sfit) )
    case 'spheremodel'
        C = reshape(sfit.Center,3,1);
        r = sfit.Radius;
    case 'double'
        if numel(sfit) == 4
            C = reshape(sfit(1:3),3,1);
            r = sfit(4);
        else
            msg = msgErr;
        end
    case 'single'
        if numel(sfit) == 4
            C = reshape(sfit(1:3),3,1);
            r = sfit(4);
        else
            msg = msgErr;
        end
    otherwise
        msg = msgErr;
end

if ~isempty(msg)
    error(msg);
end

%% Plot sphere
[x,y,z] = sphere(N);
[m,n] = size(x);

X_s = [reshape(x,1,[]); reshape(y,1,[]); reshape(z,1,[])];
% Apply radius
X_s = r*X_s;
% Make points homogeneous
X_s(4,:) = 1;

% Define transformation
H_s2a = eye(4);
H_s2a(1:3,4) = C;

if nargout < 2
    X_a = H_s2a*X_s;
    x = reshape(X_a(1,:),m,n);
    y = reshape(X_a(2,:),m,n);
    z = reshape(X_a(3,:),m,n);
    
    %srf = surf(axs,x,y,z,'Visible','off');
    srf = surf(x,y,z,'Visible','off','Parent',axs);
    ptc = surf2patch(srf);
    delete(srf);
    %p_a = patch(axs,ptc,'EdgeColor','None','FaceColor','b','FaceAlpha',0.5);
    p_a = patch(ptc,'EdgeColor','None','FaceColor','b','FaceAlpha',0.5,...
        'Parent',axs);
    if nargout == 1
        varargout{1} = p_a;
    end
else
    x = reshape(X_s(1,:),m,n);
    y = reshape(X_s(2,:),m,n);
    z = reshape(X_s(3,:),m,n);
    %srf = surf(axs,x,y,z,'Visible','off');
    srf = surf(x,y,z,'Visible','off','Parent',axs);
    ptc = surf2patch(srf);
    delete(srf);
    h_s2a = triad('Parent',axs,'Scale',r,'LineWidth',1.5,'Matrix',H_s2a);
    p_s = patch(ptc,'EdgeColor','None','FaceColor','b','FaceAlpha',0.5,...
        'Parent',h_s2a);
    varargout{1} = p_s;
    varargout{2} = h_s2a;
end