%% Description
% Compute the signal norm

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function nrm = signal_norm(u,dom)
% u         (scalar) signal
% dom       domain: 1 - time, 2 - freq

if dom == 1 
    % Time domain: 
    %   nrm = sqrt(\sum_k(u_k^2)/N)
    [a,b] = size(u);
    if a>b
        nrm = sqrt(u'*u/length(u));
    else
        nrm = sqrt(u*u'/length(u));
    end
else
    % Freq domain: 
    % 	nrm = sqrt(\sum_k |U[k]|^2 for all k)
    %           or 
    %   nrm = sqrt(|U[0]|^2+ 2\sum_{k=1}^N/2 |U[k]|^2)
    % The latter is used here
    [a,b] = size(u);
    if a > b
        nrm = sqrt(u(1)^2+2*u(2:round(a/2))'*u(2:round(a/2)));
    else
        nrm = sqrt(u(1)^2+2*u(2:round(b/2))*u(2:round(b/2))');
    end
end