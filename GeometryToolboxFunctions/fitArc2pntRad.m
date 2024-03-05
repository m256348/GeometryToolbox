function afit = fitArc2pntRad(X,r)
% FITARC2PNTRAD fits a circular arc to two, two dimensional points and a
% raduis.
%
%   afit = fitArc2pntRad(X,r)
%
%   Input(s)
%       X - 2x2 array containing points
%           pnt(:,1) - x/y coordinate of point 1
%           pnt(:,2) - x/y coordinate of point 2
%       r - scalar, non-zero value defining radius
%
%   Output(s)
%       afit - structured array containing the following fields
%           afit.Center    - 3x1 center of the arc
%           afit.Rotation  - 3x3 orientation of the arc (rotation is 
%                              about the z-direction)
%           afit.Radius    - radius of the arc
%           afit.AngleLims - Nx2 array containing the bounds of the 
%                            angles used to define the arc.
%                AngleLims(i,:) - lower and upper bounds of the angle
%                                 defining the ith arc intersection.
%
%   NOTE:
%       This function returns 2D arcs meaning the center of the arc will
%       have a z-component of 0.
%
%   M. Kutzer, 04Mar2024, USNA

debugOn = true;

%% Check Input(s)
narginchk(2,2)

% Check points
[m,n] = size(X);
if m ~= 2 || n ~= 2
    error('Specified points must be provided as a 2x2 array.');
end

% Check radius
if numel(r) ~= 1 || ~isreal(r) || r == 0
    error('Radius must be a scalar, non-zero value.');
end


%% Fit arc
% Define vector between points
v = diff(X,1,2);    % Vector pointing from pnt1 to pnt2
v_hat = v./norm(v); % Unit vector

% Check for ill-defined radius
if abs(r) < norm(v)/2
    error('Radius must be at least 1/2 the distance between the points.');
end

% Define normal vector
n = nCross(v_hat);      % Normal vector
n_hat = n./norm(n); % Unit vector

% Define midpoint between points
X_m = mean(X,2);    % Midpoint

% Define center of circle
d = sign(r)*sqrt(r^2 - (norm(v)/2)^2);
X_c = X_m + d*n_hat;

% Define angles
phi(1) = atan2(X(2,1) - X_c(2), X(1,1) - X_c(1));
a = X(:,1) - X_c;
b = X(:,2) - X_c;
dphi = -sign(r)*acos( dot(a,b)/(norm(a)*norm(b)) );
phi(2) = phi(1) + dphi;

%% Package output(s)
afit.Center    = [X_c; 0];
afit.Rotation  = eye(3);
afit.Radius    = abs(r);
afit.AngleLims = phi;

%% Debug plot
if debugOn
    fig = figure('Name',mfilename);
    axs = axes('Parent',fig,'NextPlot','add','DataAspectRatio',[1 1 1]);
    
    pnt = plot(axs,X(1,:),X(2,:),'ob');
    cnt = plot(axs,afit.Center(1,:),afit.Center(2,:),'xg');
    
    for i = 1:2
        con(i) = plot(axs,...
            [afit.Center(1),X(1,i)],...
            [afit.Center(2),X(2,i)],':k');
        txt(i) = text(axs,X(1,i),X(2,i),sprintf('p_{%d}',i),...
            'HorizontalAlignment','Right','VerticalAlignment','Top');
    end
    
    XX = interpArc(afit,1000);
    xx = XX(1,:);
    yy = XX(2,:);
    
    crc = plot(axs,xx,yy,'b');
    crc0 = plot(axs,xx(1),yy(1),'*b');
    crc1 = plot(axs,xx(end),yy(end),'xb');
end
