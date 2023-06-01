%% Description
% Example of Lemma 18 in the paper [] on model order
% reduction by moment matching for convergent Lur'e-type models.

% The hinfnorm of the reduced-order model sysr should not exceed that of
% the full-order model sys. The system sysr_rand has been reduced with a
% random matrix G, so there is no guarantee that its hinfnorm remains below
% that of the full-order model sys. 

% Author: Fahim Shakib
% Date:   February. 21, 2022
% Email:  m.f.shakib@tue.nl

%% Initialization
clear all; clc

% Maximum state dimension full-order model
imax    = 10;

%% Take random full-order model
n       = randi(imax);
sys0    = rss(n);

% Compute its Hinf norm - gamma is the infinity norm
gamma = norm(sys0,inf);

% Find the matrix Pbar
Aplus = sys0.A+1/gamma*sys0.B*sys0.C;
Amin  = sys0.A-1/gamma*sys0.B*sys0.C;

Pbar = sdpvar(n);
LMI = Pbar>=eye(n)*eps;
LMI = [LMI, Pbar*Aplus+Aplus'*Pbar<=-eye(n)*eps];
LMI = [LMI, Pbar*Amin+Amin'*Pbar<=-eye(n)*eps];

sol = optimize(LMI);

if ~double(any(~checkset(LMI)>0))
    Pbar = double(Pbar);
else
    display('LMIs infeasible')
    return
end

%% Take random interpolation points
flg = 1;
while flg
    v = randi(imax);
    flg = ~(v < n);
end

S = randn(v);

% Select L such that (S,L) is observable
flg = 1;
while flg
    L   = randn(1,v);
    flg = ~(rank(obsv(S,L)) == v);
end

%% Compute Pi
Pi = lyap(sys0.A,-S,sys0.B*L);

%% Compute G according to Lemma 18
G = (Pi'*Pbar*Pi)\Pi'*Pbar*sys0.B;
Grand = randn(size(G));

%% Construct reduced order model
F       = S-G*L;
Frand   = S-Grand*L;
H       = sys0.C*Pi;

sysr        = ss(F,G,H,0);
sysr_rand   = ss(Frand,Grand,H,0);

% Compute its Hinf norm
gamma_red       = norm(sysr,inf);
gamma_red_rand  = norm(sysr_rand,inf);

%% Print results
display(['Hinf full-order model ' num2str(gamma)])
display(['Hinf reduced-order model ' num2str(gamma_red)])
display(['Hinf random reduced-order model ' num2str(gamma_red_rand)])

%% Show results in Bode plot
figure
bode(sys0)
hold all
bode(sysr)
bode(sysr_rand)
legend('Full-order model','Reduced-order model','Random reduced-order model')
