%% Description
% Modified cost function with a penalty function to handle MIMO models

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function [J,sysr] = cost_fnc_FRF_penalty_MIMO(G,H,S,L,mag0,wout0,v,gamma)

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

% Define reduced-order LTI model
sysr     = ss(S-G*L,G,H,0);

% Compute Hinf(sysred)*gamma. Alternatively, this could also be computed
% via the LMIs in Theorem 19. However, this implementation is faster
cond     = norm(sysr(1,2),inf)*gamma;
delta    = 1e2;
% Assign a value for c
if cond<=1
    c    = 1;
else
    c    = 1+(cond-1)/delta;
end

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
% Multiply by c for a penalty when cond > 1 (,i.e., reduced-order model
% does not satisfy Hinf(sys_red)< 1/gamma)
J = J*c;