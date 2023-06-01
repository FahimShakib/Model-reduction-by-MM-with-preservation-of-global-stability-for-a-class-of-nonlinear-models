%% Description
% Simulate the steady-state response of wiener type of nonlinear models

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function z = sim_Wiener(u,t,sys,varphi)
% inputs
%   u       input signal
%   t       Time vector
%   sys     LTI part of Wiener dynamics
%   varphi  Nonlinearity
%
% Outputs
%   z       Steady-state model output


% Simulate for P periods, only keep the last one as the steady-state
P = 4;

% Extend time and input sequence to P periods
ts      = t(2)-t(1);
T       = t(end)+ts;
t_ext   = 0:ts:T*P-ts;   % Time vector
u_ext   = [];
for i = 1:P
    u_ext = [u_ext u];
end

[~,~,x_ext] = lsim(sys,u_ext,t_ext);
x           = x_ext(end-length(u)+1:end,:);
z           = varphi(x);