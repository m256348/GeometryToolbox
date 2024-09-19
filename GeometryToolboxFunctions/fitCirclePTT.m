function [cfit,Xint_cfit] = fitCirclePTT(X,abc1,abc2,varargin)
% FITCIRCLEPTT fits a circle to a point and two lines assuming lines are
% tangent and coincident to the circle.
%   [cfit,Xint_cfit] = fitCirclePTT(X,abc1,abc2)
%   ___ = fitCirclePTT(___,debug)
%   ___ = fitCirclePTT(___,ZERO)
%
%   Input(s)
%       X - 2x1 array defining 2D point on the circle
%    abc1 - 1x3 array defining the coefficients to a line that is both
%           coincident and tangent to the circle
%    abc2 - 1x3 array defining the coefficients to a line that is both
%           coincident and tangent to the circle
%   debug - [OPTIONAL] logical scalar indicating whether or not to show 
%           debug plots. Default value is false.
%    ZERO - [OPTIONAL] positive value that is sufficiently close to zero 
%           to be assumed zero (e.g. ZERO = 1e-8). If ZERO is not 
%           specified or if ZERO = [], a default value is used.
%
%   Output(s)
%        cfit - 2-element structured array containing the following fields
%           cfit.Center - 3x1 center of the circle (z-coordinate is 0)
%           cfit.Normal - 3x1 normal to the circle ([0; 0; 1])
%           cfit.Radius - radius of the circle
%   Xint_cfit - 2-element cell array defining the points of intersection
%               between the lines and the circle
%       Xint_cfit{i}(:,j) - x/y intersection of circle i with line j
%
%   NOTE: cfit fields are to retain compatibility with fitCircle.m
%
%   See also fitCircle fitPlane
%
%   M. Kutzer, 14May2024, USNA

% Update(s)
%   04Sep2024 - Updated debug plots, special case considerations, and
%               assignin use for non-real parts
%   05Sep2024 - Updated for the special case where Line i (i \in {1,2}) 
%               contains the point X, two candidate circles are possible. 
%               Each of these points results from the two possible +/- 
%               orientations of Line i.
%   05Sep2024 - Updated to allow user to specify debug and ZERO.

%% Check input(s)
narginchk(3,5);

if numel(X) ~= 2 || ~isnumeric(X)
    error('Point on circle must be 2D.');
end

if numel(abc1) ~= 3 || ~isnumeric(abc1)
    error('Line 1 must be defined using three constants.');
end

if numel(abc2) ~= 3 || ~isnumeric(abc1)
    error('Line 1 must be defined using three constants.');
end

%% Reshape input(s)
X = reshape(X,[],1);
abc1 = reshape(abc1,1,[]);
abc2 = reshape(abc2,1,[]);

% Check if lines are the same
if rank( [abc1; abc2] ) < 2
    error('Tangent lines are the same.')
end

% Set defaults
debug = false;
ZERO = 1e-8;
if nargin > 3
    % Parse values
    for i = 1:numel(varargin)
        if numel(varargin{i}) > 1
            error('Optional values must be scalars.');
        end
        if numel(varargin{i}) == 0
            continue
        end

        if islogical(varargin{i})
            debug = varargin{i};
        else
            ZERO = varargin{i};
        end
    end
end

if ZERO < 0
    error('ZERO must be a positive value.');
end

%% Initialize output
cfit.Center = [];
cfit.Normal = [0;0;1];
cfit.Radius = [];

%% Check for special case
tfSpecialCase = false(1,2);

% Point lies on Line 1
if abs(abc1*[X; 1]) < ZERO
    tfSpecialCase(1) = true;

    % Consider both line orientations
    abc1 = ...
        [abc1;
        -abc1];
end

% Point lies on Line 2
if abs(abc2*[X; 1]) < ZERO
    tfSpecialCase(2) = true;

    % Consider both line orientations
    abc2 = ...
        [abc2;
        -abc2];
end

% Point lies on both Line 1 and Line 2
if all(tfSpecialCase,'all')
    msg = sprintf([...
        'Tangent lines intersect at the point given.\n',...
        'The resultant circle has a radius of 0.']);
    warning(msg);
    cfit.Center = [X; 0];
    cfit.Radius = 0;
    Xint_cfit{1} = X;
    return
end


%% Create debug plot
if debug
    % Create figure and axes
    fig = figure('Name','fitCirclePTT.m, debug = true');
    axs = axes('Parent',fig,'DataAspectRatio',[1 1 1],'NextPlot','add');

    % Plot point
    plt_X = plot(axs,X(1),X(2),'ok','MarkerFaceColor','k');
end

%% Orient lines
% Place point in positive half space
if ~tfSpecialCase(1)
    if abc1*[X; 1] < 0
        abc1 = -abc1;
    end
end

if ~tfSpecialCase(2)
    if abc2*[X; 1] < 0
        abc2 = -abc2;
    end
end

%% Calculate intersect between lines
% -> Point of intersection is the same regardless of line orientation
Xint = intersectLineLine(abc1(1,:),abc2(1,:));

% Define distance between points
dX = norm(X-Xint);
xx = Xint(1) + dX*[-1,1];
yy = Xint(2) + dX*[-1,1];

if debug
    % Plot intersection
    plt_Xint = plot(axs,Xint(1),Xint(2),'ok');
    
    % Plot xx & yy
    for i = 1:numel(xx)
        for j = 1:numel(yy)
            plt_xxyy(i,j) = plot(axs,xx(i),yy(j),'oc');
            txt_xxyy(i,j) = text(axs,xx(i),yy(j),sprintf('[%d,%d]',i,j),...
                'HorizontalAlignment','center','VerticalAlignment','Top');
        end
    end
end

%% Calculate angle between lines
X1 = line2points(abc1(1,:),xx,yy,ZERO);
X2 = line2points(abc2(1,:),xx,yy,ZERO);

if debug
    lsymbols = '+x';
    for i = 1:2
        plot(axs,X1(1,i),X1(2,i),['r',lsymbols(i)]);
        plot(axs,X2(1,i),X2(2,i),['g',lsymbols(i)]);
    end
end

% Keep points in positive half space
% -> Point on Line 1
for i = 1:size(abc2,1)
    tf = abc2(i,:)*[X1; 1 1] > 0;
    X1_TMP(:,i) = X1(:,tf);
end
% -> Point on Line 2
for i = 1:size(abc1,1)
    tf = abc1(i,:)*[X2; 1 1] > 0;
    X2_TMP(:,i) = X2(:,tf);
end
X1 = X1_TMP;
X2 = X2_TMP;

if debug
    plot(axs,X1(1,:),X1(2,:),'ro');
    plot(axs,X2(1,:),X2(2,:),'go');
end

%% Calculate angles
% -> Angle to Line 1
phi1 = atan2(X1(2,:)-Xint(2),X1(1,:)-Xint(1));
% -> Angle to Line 2
phi2 = atan2(X2(2,:)-Xint(2),X2(1,:)-Xint(1));

% Account for special case
if tfSpecialCase(1)
    phi1 = repmat(phi1,1,2);
end
if tfSpecialCase(2)
    phi2 = repmat(phi2,1,2);
end

% Wrap angle(s) to ensure an arc of less than pi
tf = abs(phi2 - phi1) >= pi;
phi1(tf) = wrapTo2Pi(phi1(tf));
phi2(tf) = wrapTo2Pi(phi2(tf));

% Calculate arc angle 
dphi = abs( wrapToPi(phi2 - phi1) );
phiB = phi1 + sign(phi2-phi1).*dphi./2;

% fprintf('phi1 = %.2f, phi2 = %.2f, phiB = %.2f\n',...
%     rad2deg(phi1),rad2deg(phi2),rad2deg(phiB));

if debug
    for j = 1:numel(phiB)
        afit(1,j).Center    = [Xint;0];
        afit(1,j).Rotation  = eye(3);
        afit(1,j).Radius    = 0.8*norm(X-Xint)/2;
        afit(1,j).AngleLims = [phi1(j),phi2(j)];

        afit(2,j).Center    = [Xint;0];
        afit(2,j).Rotation  = eye(3);
        afit(2,j).Radius    = norm(X-Xint)/2;
        afit(2,j).AngleLims = [phi1(j),phiB(j)];

        afit(3,j).Center    = [Xint;0];
        afit(3,j).Rotation  = eye(3);
        afit(3,j).Radius    = norm(X-Xint)/2;
        afit(3,j).AngleLims = [phi2(j),phiB(j)];
    end
    
    colors = 'krg';
    lwidth = [...
        0.5,1.0,0.5;...
        1.5,2.0,1.5];
    for j = 1:size(afit,2)
        for i = 1:size(afit,1)
            pltArc(i,j) = plotArc(axs,afit(i,j));
            set(pltArc(i,j),'Color',colors(i),'LineWidth',lwidth(j,i));
        end
    end

    % Define "line B"
    markers = {'--',':'};
    for i = 1:numel(phiB)
        abcB(i,:) = fitPlane(Xint + [0,dX*cos(phiB(i)); 0,dX*sin(phiB(i))]);
        plt_abcB(i) = plotLine(axs,abcB(i,:),[markers{i},'m'],...
            'LineWidth',1);
    end
end

%% Display lines
if debug
    % Plot line 1
    plt_abc1 = plotLine(axs,abc1(1,:),'-r','LineWidth',1);

    % Plot line 2
    plt_abc2 = plotLine(axs,abc2(1,:),'-g','LineWidth',1);

    % Force axis limits to "tight"
    axis(axs,'tight');
end

%% Solve quadratic equation
h = [];
phi1_TMP = [];
phi2_TMP = [];
dphi_TMP = [];
phiB_TMP = [];
for i = 1:numel(phiB)
    % r^2 = (x - (h*cos(phiB)+x_int))^2 + (y - (h*sin(phiB)+y_int)^2)
    % r^2 = h^2 * sin(dphi/2)
    x = X(1);
    y = X(2);
    x_int = Xint(1);
    y_int = Xint(2);
    cpB = cos(phiB(i));
    spB = sin(phiB(i));
    sp  = sin(dphi(i)/2);
    % solve((x - (h*cpB+x_int))^2 + (y - (h*spB+y_int))^2 == (h*sp)^2,h)

    h_TMP(1) = ...
        ((- cpB^2*y^2 + 2*cpB^2*y*y_int - cpB^2*y_int^2 + 2*cpB*spB*x*y ...
        - 2*cpB*spB*x*y_int - 2*cpB*spB*x_int*y + 2*cpB*spB*x_int*y_int ...
        + sp^2*x^2 - 2*sp^2*x*x_int + sp^2*x_int^2 + sp^2*y^2 ...
        - 2*sp^2*y*y_int + sp^2*y_int^2 - spB^2*x^2 + 2*spB^2*x*x_int ...
        - spB^2*x_int^2)^(1/2) + cpB*x - cpB*x_int + spB*y - spB*y_int)...
        /(cpB^2 - sp^2 + spB^2);
    h_TMP(2) = ...
        -((- cpB^2*y^2 + 2*cpB^2*y*y_int - cpB^2*y_int^2 + 2*cpB*spB*x*y ...
        - 2*cpB*spB*x*y_int - 2*cpB*spB*x_int*y + 2*cpB*spB*x_int*y_int ...
        + sp^2*x^2 - 2*sp^2*x*x_int + sp^2*x_int^2 + sp^2*y^2 ...
        - 2*sp^2*y*y_int + sp^2*y_int^2 - spB^2*x^2 + 2*spB^2*x*x_int ...
        - spB^2*x_int^2)^(1/2) - cpB*x + cpB*x_int - spB*y + spB*y_int)...
        /(cpB^2 - sp^2 + spB^2);

    % Check for special case (non-real values)
    % TODO - validate this approach
    if any(tfSpecialCase)
        %if abs(diff( real(h_TMP) )) < ZERO
        h_TMP = sqrt( prod(h_TMP) );
        %end
    end
    
    % Combine results
    h = [h,h_TMP];
    phi1_TMP = [phi1_TMP,repmat(phi1(i),1,numel(h_TMP))];
    phi2_TMP = [phi2_TMP,repmat(phi2(i),1,numel(h_TMP))];
    dphi_TMP = [dphi_TMP,repmat(dphi(i),1,numel(h_TMP))];
    phiB_TMP = [phiB_TMP,repmat(phiB(i),1,numel(h_TMP))];
end

% Update angles to angles matched with h
phi1 = phi1_TMP;
phi2 = phi2_TMP;
dphi = dphi_TMP;
phiB = phiB_TMP;

%% Recover circle parameters
for i = 1:numel(h)
    X0(1,i) = h(i)*cos(phiB(i)) + Xint(1);
    X0(2,i) = h(i)*sin(phiB(i)) + Xint(2);

    imagVals = imag(X0(:,i));
    if max(abs(imagVals)) > ZERO
        % Center has imaginary components
        % -> Define message to user
        str = sprintf([...
            'cfit(%d).Center has an imaginary part:\n',...
            '\tcfit(%d).Center = [%.8f,%.8f] + [%.8f,%.8f]i\n',...
            '\t -> Forcing real values.\n',...
            '[cfitTST,Xint_cfitTST] = fitCirclePTT(debug_fitCirclePTT(%d).X,...\n',...
            '\tdebug_fitCirclePTT(%d).abc1,debug_fitCirclePTT(%d).abc2,true);\n'],...
            i,i,real(X0(:,i)),imagVals);
            
        % -> Package information for debugging
        debug_fitCirclePTT(i).X     = X;
        debug_fitCirclePTT(i).abc1  = abc1;
        debug_fitCirclePTT(i).abc2  = abc2;
        debug_fitCirclePTT(i).x     = x;
        debug_fitCirclePTT(i).y     = y;
        debug_fitCirclePTT(i).x_int = x_int;
        debug_fitCirclePTT(i).y_int = y_int;
        %debug_fitCirclePTT(i).cpB   = cpB;
        %debug_fitCirclePTT(i).spB   = spB;
        %debug_fitCirclePTT(i).sp    = sp;
        debug_fitCirclePTT(i).phiB  = phiB(i);
        debug_fitCirclePTT(i).dphi  = dphi(i);
        debug_fitCirclePTT(i).h     = h(i);
        debug_fitCirclePTT(i).X0    = X0(:,i);
        % -> Assign debug values to base workspace
        assignin('base','debug_fitCirclePTT',debug_fitCirclePTT);
        warning(str);
    end
    X0(:,i) = real(X0(:,i));

    r(i) = norm(X - X0(:,i));
    d(i) = sqrt(h(i)^2 - r(i)^2);

    cfit(i).Center = [X0(:,i); 0];
    cfit(i).Normal = [0; 0; 1];
    cfit(i).Radius = r(i);
    
    Xint_cfit{i}(1,1) = real( d(i)*cos(phi1(i)) + Xint(1) );
    Xint_cfit{i}(2,1) = real( d(i)*sin(phi1(i)) + Xint(2) );

    Xint_cfit{i}(1,2) = real( d(i)*cos(phi2(i)) + Xint(1) );
    Xint_cfit{i}(2,2) = real( d(i)*sin(phi2(i)) + Xint(2) );

    if debug
        % Define markers
        markers = '+x';
        % Define colors
        colors = 'bc';
        % Define linewidths
        linewidths = [3,2];
        % Plot center point
        pltX0(i) = plot(axs,X0(1,i),X0(2,i),[markers(i),colors(i)]);
        
        % Plot arcs
        ptcCrc(i) = plotCircle(axs,cfit(i));
        set(ptcCrc(i),'LineWidth',0.5,'FaceColor',colors(i),...
            'EdgeColor',colors(i),'FaceAlpha',0.1);
        
        % Plot h
        plt_h(i) = plot(axs,[Xint(1),X0(1,i)],[Xint(2),X0(2,i)],colors(i),...
            'LineWidth',linewidths(i));

        % Update Line plots
        delete(plt_abc1);
        delete(plt_abc2);
        delete(plt_abcB);
        % -> Plot line 1
        plt_abc1 = plotLine(axs,abc1(1,:),'-r','LineWidth',1);
        % -> Plot line 2
        plt_abc2 = plotLine(axs,abc2(1,:),'-g','LineWidth',1);
        % -> Plot line B
        markers = {'--',':'};
        for j = 1:size(abcB,1)
            plt_abcB(j) = plotLine(axs,abcB(j,:),[markers{j},'m'],...
                'LineWidth',1);
        end
    end

end

if debug
    % Create legend
    if numel(plt_h) == 1
        lgnd = legend(...
            [plt_X,plt_abc1,plt_abc2,pltArc(1,1),pltArc(2,1),...
            pltArc(3,1),plt_abcB,pltX0(1),ptcCrc(1),plt_h(1)],...
            {'Point X','Line 1','Line 2','\phi_1 to \phi_2','\phi_1 to \phi_B',...
            '\phi_2 to \phi_B','Line B','Center 1 (X_0)','Circle 1','h_1'},...
            'Location','northeastoutside');
    else
        lgnd = legend(...
            [plt_X,plt_abc1,plt_abc2,pltArc(1,1),pltArc(2,1),pltArc(3,1),...
            plt_abcB(1),pltX0(1),ptcCrc(1),plt_h(1),pltX0(2),ptcCrc(2),plt_h(2)],...
            {'Point X','Line 1','Line 2','\phi_1 to \phi_2','\phi_1 to \phi_B',...
            '\phi_2 to \phi_B','Line B','Center 1 (X_0)','Circle 1','h_1',...
            'Center 2 (X_0)','Circle 2','h_2'},'Location','northeastoutside');
    end
end
