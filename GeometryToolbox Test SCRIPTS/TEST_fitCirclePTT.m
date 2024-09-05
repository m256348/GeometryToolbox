%% TEST_fitCirclePTT
% Test fitCirclePTT
%
%   M. Kutzer, 14May2024, USNA
clear all
close all
clc

%% Create random data set
%load('TestCase_fitCirclePTT.mat');

% Generate random points
pnts = 10*(rand(2,5)-0.5);
% Select first point
desc{1} = 'Free Point';
X = pnts(:,1);
% Fit planes
abc(1,:) = fitPlane(pnts(:,2:3));
abc(2,:) = fitPlane(pnts(:,4:5));

% Define special cases
desc{2} = 'Point on Line 1';
X(:,2) = proj2line(abc(1,:),X(:,1));
desc{3} = 'Point on Line 2';
X(:,3) = proj2line(abc(2,:),X(:,1));
desc{4} = 'Point on Line 1 & 2';
X(:,4) = intersectLineLine(abc(1,:),abc(2,:));

%% Test functions
ZERO = 1e-8;
for i = 1:numel(desc)
    fprintf('---------- %s ----------\n',desc{i});
    % Run function
    [cfit,Xint_cfit] = fitCirclePTT(X(:,i),abc(1,:),abc(2,:),true);
    % Test results
    for j = 1:numel(cfit)
        % Check if intersections lie on circle
        Xint = Xint_cfit{j};
        val = evalCfit(cfit(j),Xint);
        tf = abs(val) > ZERO;
        disp(tf);
        [a,b] = find(tf);
        for k = 1:numel(a)
            fprintf('Fit %d, Point %d - ',j,a(k))
            switch b(k)
                case 1
                    fprintf('Not on circle.\n');
                case 2
                    fprintf('Not on the plane of the circle\n');
            end
        end
        
        % Check if intersections lie on line
        Xint(3,:) = 1;
        val = abc*Xint;
        tf = abs(val) > ZERO;
        disp(tf);
        [a,b] = find(tf);
        for k = 1:numel(a)
            fprintf('Line %d, Point %d - Is not on line\n',a(k),b(k))
        end
    end
end


