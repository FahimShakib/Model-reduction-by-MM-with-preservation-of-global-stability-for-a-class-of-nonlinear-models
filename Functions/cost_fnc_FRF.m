%% Description
% Modified cost function with a penalty function

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function [J,sysr] = cost_fnc_FRF(theta,H,S,L,mag0,wout0,v,gamma)

% Inputs:
%   theta     To-be-optimized parameters
%   H         Output matrix of reduced-order model
%   S         Signal generator system matrix
%   L         Signal generator output matrix
%   mag0      Magnitude of Sigma
%   wout0     Frequencies corresponding to mag0
%   v         Order of reduced-order LTI models
%   gamma     Maximum slope of nonlinearity

% Outputs:
%   J         Value of modified cost function
%   sysr      Reduced-order model with updated matrices

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

% Define reduced-order LTI model
sysr     = ss(FROM,[G1ROM G2ROM],[H1ROM;H2ROM],0);

% Compute Hinf(sysred)*gamma. Alternatively, this could also be computed
% via the LMIs in Theorem 19. However, this implementation is faster
cond     = norm(sysr(1,2),inf)*gamma;

% Compute FRF of sysr at wout0
[magr,~,~] = bode(sysr,wout0);

% Compute J
e = mag0-magr;
J = 0;
for i = 1:2
    for k = 1:2
        etmp = squeeze(e(i,k,:));
        J = J + etmp'*etmp/length(etmp);
    end
end