function varargout = GeometryToolboxVer
% GEOMETRYTOOLBOXVER displays the Piecewise Polynomial Toolbox information.
%   GEOMETRYTOOLBOXVER displays the information to the command prompt.
%
%   A = GEOMETRYTOOLBOXVER returns in A the sorted struct array of  
%   version information for the Piecewise Polynomial Toolbox.
%     The definition of struct A is:
%             A.Name      : toolbox name
%             A.Version   : toolbox version number
%             A.Release   : toolbox release string
%             A.Date      : toolbox release date
%
%   M. Kutzer 28Aug2018, USNA

A.Name = 'Geometry Toolbox';
A.Version = '1.0.0';
A.Release = '(R2017b)';
A.Date = '28-Aug-2018';
A.URLVer = 1;

msg{1} = sprintf('MATLAB %s Version: %s %s',A.Name, A.Version, A.Release);
msg{2} = sprintf('Release Date: %s',A.Date);

n = 0;
for i = 1:numel(msg)
    n = max( [n,numel(msg{i})] );
end

fprintf('%s\n',repmat('-',1,n));
for i = 1:numel(msg)
    fprintf('%s\n',msg{i});
end
fprintf('%s\n',repmat('-',1,n));

if nargout == 1
    varargout{1} = A;
end