%% Description
% One-sided spring nonlinearity

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function [u] = onesidedspring(y,gamma)
% inputs
%   y       input signal
%   gamma   slope of spring
%
% Outputs
%   u       output signal


u = max(gamma*y,0);
