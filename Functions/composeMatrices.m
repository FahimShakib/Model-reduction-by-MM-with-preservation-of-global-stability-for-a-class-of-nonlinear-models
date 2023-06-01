%% Description
% Composes the sparse stiffness and mass matrices based on the given 
% number of elements.

% Author: Fahim Shakib (based on files received from Luuk Poort, M.Sc.)
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

function [K,M] = composeMatrices(Ne, L, E, I, rho, A)
    Nn  = Ne+1;           % Number of nodes
    Nd  = 2*Nn;           % Number of DoF
    l   = L/Ne;           % Beam element length
    K   = sparse(Nd,Nd);  % Stiffnes matrix
    M   = sparse(Nd,Nd);  % Mass matrix

    % Elemental stiffness and mass matrices
    Ke = E*I/l^3 * [
        12,     6*l,    -12,    6*l;
        6*l,    4*l^2,  -6*l,   2*l^2;
        -12,    -6*l,   12,     -6*l;
        6*l,    2*l^2,  -6*l,   4*l^2];
    Me = rho*A*l/420 * [
        156,    22*l,   54,     -13*l;
        22*l,   4*l^2,  13*l,   -3*l^2;
        54,     13*l,   156,    -22*l;
        -13*l,  -3*l^2, -22*l,  4*l^2];

    % Compose the mass and stiffness matrices
    for ii = 1:Ne
        indx = ii*2-1:ii*2+2;
        K(indx,indx) = K(indx,indx) + Ke;
        M(indx,indx) = M(indx,indx) + Me;
    end
end