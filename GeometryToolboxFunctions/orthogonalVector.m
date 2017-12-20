function v = orthogonalVector(u)
% ORTHOGONALVECTOR finds a unit vector orthogonal to the vector provided.
%   v = ORTHOGONALVECTOR(u) finds a unit vector orthogonal to the vector
%    provided
%
%   See also
%
%   M. Kutzer, 20Dec2017, USNA

%% Check Inputs
narginchk(1,1);

[m,n] = size(u);
if m ~= 3 || n ~= 1
    error('Vector must be a 3x1.');
end

%% Find orthogonal vector 
so = soBasis(m);
for i = 1:numel(so)
    v(:,i) = so{i}*u;
end

v_mag = sqrt( sum(v.^2,1) );
[~,idx] = sort(v_mag,'descend');

v = v(:,idx(1))./v_mag(idx(1));