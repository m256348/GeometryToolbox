%% SCRIPT_Solve_fitCirclePTT
% !!! This is not a useful solution !!!
% -> Please see fitCirclePTT

clear all
close all
clc

%% Define symbolic terms
syms a1 b1 c1 a2 b2 c2 x y

X = [x;y];
abc1 = [a1,b1,c1];
abc2 = [a2,b2,c2];

%abP1 = nCross(abc1(1:2));
abP1 = [abc1(2),-abc1(1)];
%abP2 = nCross(abc2(1:2));
abP2 = [abc2(2),-abc2(1)];

aP1 = abP1(1);
bP1 = abP1(2);

aP2 = abP2(1);
bP2 = abP2(2);

%% Define unknown symbolic terms
syms cP1 cP2

X0 = -([aP1, bP1; aP2, bP2])^(-1) * [cP1; cP2];
X1 = -([a1, b1; aP1, bP1])^(-1) * [c1; cP1];
X2 = -([a2, b2; aP2, bP2])^(-1) * [c2; cP2];

[cP1_out,params,conds] = solve(...
    norm(X - X0) == norm(X1 - X0),...
    cP1,...
    'ReturnConditions',true );
conds 

return
%%
X0_out = subs(X0,cP1,cP1_out);
X2_out = subs(X2,cP1,cP1_out);

[cP2_out,params,conds] = solve(...
    norm(X - X0_out) == norm(X2_out - X0_out),...
    cp2,...
    'ReturnConditions',true );
conds 

X0_out = subs(X0,cP2,cP2_out);
X1_out = subs(X1,cP2,cP2_out);

[cP1_out,params,conds] = solve(...
    norm(X - X0) == norm(X1 - X0),...
    cP1,...
    'ReturnConditions',true );
conds