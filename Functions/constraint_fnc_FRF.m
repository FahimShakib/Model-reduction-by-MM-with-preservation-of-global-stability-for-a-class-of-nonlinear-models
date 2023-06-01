%% Description
% Constraints eqn (27) 

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function [cineq,ceq] = constraint_fnc_FRF(theta,H,S,L,v,gamma)

% Inputs:
%   theta     To-be-optimized parameters
%   H         Output matrix of reduced-order model
%   S         Signal generator system matrix
%   L         Signal generator output matrix
%   v         Order of reduced-order LTI models
%   gamma     Maximum slope of nonlinearity

% Outputs:
%   cineq     -1 if feasible, +1 if infeasible
%   ceq       [];

% Define G and F
FROM = [];
kk = 0;
for i = 1:2
    for k = 1:2
        GG{i,k} = theta(kk+1:kk+v{i,k});
        kk      = kk + v{i,k};
        F{i,k}  = S{i,k} - GG{i,k}*L{i,k};
        FROM    = blkdiag(FROM,F{i,k});
    end
end

G1ROM = [GG{1,1}; GG{1,2}*0; GG{2,1}; GG{2,2}*0];
G2ROM = [GG{1,1}*0; GG{1,2}; GG{2,1}*0; GG{2,2}];

H1ROM = [H{1,1} H{1,2} H{2,1}*0 H{2,2}*0];
H2ROM = [H{1,1}*0 H{1,2}*0 H{2,1} H{2,2}];

F_plus = FROM+gamma*G2ROM*H1ROM;
F_min  = FROM-gamma*G2ROM*H1ROM;

vtotal = size(FROM,1);
P      = sdpvar(vtotal);
LMI    = [P >= eps*eye(vtotal)];
LMI    = [LMI, P*F_plus + F_plus'*P <= -eps*eye(vtotal)];
LMI    = [LMI, P*F_min  + F_min'*P  <= -eps*eye(vtotal)];

ops = sdpsettings;
ops.verbose = 0;
d      = solvesdp(LMI,[],ops);

feas = sum(checkset(LMI)>0) == 3;
if feas; cineq = -1; else cineq = 1; end
% If feasible, you get the output -1
% If not feasibele, you get the output one

ceq = [];