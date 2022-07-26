%% Description
% This files contains the example in the paper [] on model order
% reduction by moment matching for convergent Lur'e-type models.

% Author: Fahim Shakib
% Date:   Feb. 22, 2022
% Email:  m.f.shakib@tue.nl

% Matlab Version 9.3.0.713579 (R2017b)
% Depending software:
%   yalmip Matlab toolbox Version 20180209 (http://yalmip.github.io)
%   Mosek Version 9.1.12 (https://www.mosek.com/)

%% Initialization
clear all; clc
addpath('Functions')

%% Beam parameters
Le  = 2;            % Beam length [m] 
h   = 3e-2;         % Beam heigth [m]
w   = 5e-2;         % Beam width [m]

E   = 200e9;        % Young modulus [Pa]
rho = 7746;         % Density [kg/m^3]

A   = w*h;          % Cross section [m^2]
I   = 1/12*w^3*h;   % Second moment of inertia [m^4]

Ne  = 40;           % Number of beam elements

%% Collect beam matrices 
[K,M] = composeMatrices(Ne, Le, E, I, rho, A);
% Clamp the cantilever beam at one side
K = K(3:end,3:end);
M = M(3:end,3:end);
% Solve eigenvalue problem to obtain modes  
[U,Ld] = eig(full(K),full(M));
U = real(U);
% Add 0.1 percent modal damping for stability
Bd = (U.')\(0.1*sqrt(Ld))/U; 

%% Construct state-space model
% Retrieve size of x (size of number of second order diff. eqns)
n   = size(M,1);

% Select input location
qe  = n/2;                  % External loading
qs  = round(n/3);           % Spring force
Q   = zeros(2,n);
Q(1,qe) = 1;
Q(2,qs) = 1;

Ass = [zeros(n) eye(n); -inv(M)*K -inv(M)*Bd];
Bss = [zeros(n); inv(M)];
Css = zeros(2,2*n);
Css(1,qs) = 1;
Css(2,n-1) = 1;
Dss = 0;

sys = ss(Ass,Bss*Q',Css,Dss);

%% Simulate full-order model
% Time vector definition
ts  = 0.01;
t   = 0:ts:10-ts;

% Input signal definition
freq   = 1;
u(1,:) = sin(2*pi*freq*t);
u(2,:) = sin(2*pi*freq*t)*0;

[y,t,x] = lsim(sys,u',t);

%% Plot deflection
plot_animation = 0;
if plot_animation
    l = linspace(0,Le,Ne);
    P = 100;
    ind = round(length(t)/P);
    figure
    for i = 1:ind:length(t)
        subplot(211)
        plot(l,x(i,1:2:end/2))
        title(['Time ' num2str(t(i)) ' sec'])
        ylabel('Beam deflection [m]')
        ylim([min(min(x(:,1:2:end/2))) max(max(x(:,1:2:end/2)))]*1.1)
        subplot(212)
        plot(x(i,2:2:end/2))
        ylim([min(min(x(:,2:2:end/2))) max(max(x(:,2:2:end/2)))]*1.1)
        xlabel('Spatial coordinate [m]')
        ylabel('Beam inclination [m]')
        pause(0.1)
    end
end

%% Plot deflection at sensor location
plot_deflection = 0;
if plot_deflection
    figure
    plot(t,y)
    xlabel('Time [sec]')
    ylabel('Beam end deflection [m]')
end

%% Plot bode
plot_bode = 0;
if plot_bode
    figure
    bodemag(sys)
    xlim([1 1e4])
    grid on
end

%% Define interpolation points in [Hz]
Omega0{1,1} = [0 10.2 62.6 180];
Omega0{1,2} = [0 10.2 64.1 180];
Omega0{2,1} = [0 10.2 65.7 180];
Omega0{2,2} = [0 10.2 64.1 180];

%% Define (S,L) for moment matching and initial condition for each FRF individually
for i = 1:2
    for k = 1:2
        S{i,k} = 0;
        for wl = 2:length(Omega0{i,k})
            S{i,k} = blkdiag(S{i,k},[0 1;-1 0]*2*pi*Omega0{i,k}(wl));
        end
        sigmaS{i,k} = eig(S{i,k});
        L{i,k}      = [1; kron(ones(length(Omega0{i,k})-1,1),[1;0])]';
        tau0{i,k}   = L{i,k}';
        v{i,k}      = size(S{i,k},1);
    end
end

%% Input to MTF 
% MTF algorithm is used to simulate the Lur'e-type system below
% See [Pavlov A, Hunnekens BG, Wouw N, Nijmeijer H. Steady-state performance
% optimization for nonlinear control systems of Lur’e type. Automatica. 
% 2013 Jul 1;49(7):2087-97.] for the MTF algorithm

% MTF Settings
Tmtf      = 5;            % End time
ts        = 1e-4;           % Sampling time
tmtf      = 0:ts:Tmtf-ts;   % Time vector
max_iter  = 1000;           % # iterations in MTF
tol       = 1e-6;           % Convergence criterion of MTF

% Simulate SG1 to obtain input u
[u,~,omega] = lsim(ss(S{1,1},zeros(size(S{1,1},1),1),L{1,1},0),tmtf*0,tmtf,tau0{1,1});
u       = u';

% Plot excitation signal
plot_exc_signal = 0;
if plot_exc_signal
    figure
    plot(tmtf,u)
    xlabel('Time [sec]')
    ylabel('$u$(t)')
end

%% Define nonlinearity for full-order system
nrm_sys = norm(sys(1,2),inf);
display(['Max. \gamma^\star = ' num2str(1/nrm_sys)])

gamma   = 7.3e4;                                % Spring stiffness
NLfnc   = @(y)onesidedspring(y,gamma);          % Define one-sided spring

%% Generate data NL System
% Perfrom simulation using MTF
[y_NL,z_NL]    = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);
[y_NL0,z_NL0]  = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,@(y)(y*0));

plot_z_NL = 1;
if plot_z_NL
    figure
    plot(tmtf,z_NL,tmtf,z_NL0)
    xlabel('Time [sec]')
    ylabel('$\bar z$(t)')
end

%% Perform MOR - Solve Sylvester equations to compute CPi
% by Lyap (AX + XB + C = 0)
for i = 1:2
    for k = 1:2
        Pi{i,k}  = lyap(sys.A,-S{i,k},sys.B(:,k)*L{i,k});
        CPi{i,k} = sys.C(i,:)*Pi{i,k};
        H{i,k}   = CPi{i,k};
    end
end

%% Find initial model
% Define sdp variables
LMI = [];
for i = 1:2
    for k = 1:2
        P{i,k} = sdpvar(v{i,k});
        X{i,k} = sdpvar(v{i,k},1);
        PP{i,k} = P{i,k}*S{i,k}-X{i,k}*L{i,k};
        LMI = [LMI, P{i,k} >= eps*eye(v{i,k})];
    end
end
XX{1,2} = X{1,2}*H{1,2};

% All the v's are the same in this example
vv = v{1,1}; 

Lgp = [PP{1,1},             zeros(vv),             zeros(vv), zeros(vv);
       gamma*X{1,2}*H{1,1}, PP{1,2}+gamma*XX{1,2}, zeros(vv), zeros(vv);
       zeros(vv),           zeros(vv),             PP{2,1},   zeros(vv);
       gamma*X{2,2}*H{1,1}, gamma*X{2,2}*H{1,2},   zeros(vv), PP{2,2}];
Lgm = [PP{1,1},             zeros(vv),             zeros(vv), zeros(vv);
       -gamma*X{1,2}*H{1,1},PP{1,2}-gamma*XX{1,2}, zeros(vv), zeros(vv);
       zeros(vv),           zeros(vv),             PP{2,1},   zeros(vv);
       -gamma*X{2,2}*H{1,1},-gamma*X{2,2}*H{1,2},  zeros(vv), PP{2,2}];
   
LMI = [LMI, Lgp+Lgp'<= -eye(sum(sum(cell2mat(v))))*eps, Lgm+Lgm'<=-eye(sum(sum(cell2mat(v))))*eps];

% Solve SDP feasibility problem
sol  = solvesdp(LMI, []); 

% Retrieve solution
Finit = [];
thetainit = [];
for i = 1:2
    for k = 1:2
        P{i,k}      = double(P{i,k});
        X{i,k}      = double(X{i,k});
        G{i,k}      = inv(P{i,k})*X{i,k};
        F{i,k}      = S{i,k} - G{i,k}*L{i,k};
        Finit       = blkdiag(Finit,F{i,k});
        thetainit   = [thetainit; G{i,k}];
    end
end

G1init = [G{1,1}; G{1,2}*0; G{2,1}; G{2,2}*0];
G2init = [G{1,1}*0; G{1,2}; G{2,1}*0; G{2,2}];

H1init = [H{1,1} H{1,2} H{2,1}*0 H{2,2}*0];
H2init = [H{1,1}*0 H{1,2}*0 H{2,1} H{2,2}];

% Define intial reduced-order LTI model
sysrinit = ss(Finit,[G1init G2init],[H1init;H2init],0);

display(['Verify the convergence property: gamma*|Sys_red|_infty = ' num2str(gamma*norm(sysrinit(1,2),inf)) ' < 1'])

%% Inspect FRF for match
% Display the FRF at interpolation frequencies
for i = 1:2
    for k = 1:2
        disp(['Sys0(' num2str(i) ',' num2str(k) ')-Sys_red(' num2str(i) ',' num2str(k) '): '...
            num2str(abs(freqresp(sys(i,k)-sysrinit(i,k),2*pi*Omega0{i,k})))])
    end
end

%% Solve optimization problem
% Select frequencies to minimize cost function on
wlist = unique([logspace(-2,5,1000) Omega0{1,1} Omega0{1,2} Omega0{2,1} Omega0{2,2}])*2*pi;

% Define cost function
[mag,~,wout] = bode(sys,wlist);
options      = sdpsettings('solver','mosek','verbose',0);
J            = @(theta)cost_fnc_FRF_penalty(theta,H,S,L,mag,wout,v,gamma);

clear options
% Optimization settings
options                         = optimoptions('fminunc','Display','iter');
options.MaxFunctionEvaluations  = 10000;
options.OptimalityTolerance     = 1e-16;
% Perform optimization
thetaopt                        = fminunc(J,thetainit,options);

% Construct reduced-order LTI model
[~,sysr] = J(thetaopt);

display(['Verify the convergence property: gamma*|Sys_red|_infty = ' num2str(gamma*norm(sysr(1,2),inf)) ' < 1'])


%% Plot final bode
% Compute FRF of initial reduced-order LTI model
magrinit = bode(sysrinit,wlist);
Magrinit = squeeze(20*log10(magrinit));

% Compute FRF of full-order LTI model
mag = bode(sys,wlist);
Mag = squeeze(20*log10(mag));

% Compute FRF of initial reduced-order LTI model
magr = bode(sysr,wlist);
Magr = squeeze(20*log10(magr));

figure;
kk = 1;
for i = 1:2
    for k = 1:2
        ax{i,k} = subplot(2,2,kk);
        kk = kk + 1;
        semilogx(wlist/2/pi,squeeze(Mag(i,k,:)),'LineWidth',2);hold all
        semilogx(wlist/2/pi,squeeze(Magr(i,k,:)),'LineWidth',2,'linestyle','--');
        semilogx(wlist/2/pi,squeeze(Magrinit(i,k,:)),'LineWidth',2,'linestyle',':');
        magIP = squeeze(20*log10(bode(sys(i,k),2*pi*Omega0{i,k})));
        semilogx(Omega0{i,k},magIP,'kx','MarkerSize', 12)
    end
end
linkaxes([ax{1},ax{2},ax{3},ax{4}],'xy')
xlim([min(wlist) max(wlist)]/2/pi)
xlim([0.01 100000])
ylim([-225 -70])

subplot(ax{1,1})
ylabel('Output $y,\rho$ [dB]')
title('Input $u$')
set(gca,'fontsize', 12)
subplot(ax{1,2})
title('Input $\varphi(y),\varphi(\rho)$')
set(gca,'fontsize', 12)
subplot(ax{2,1})
ylabel('Output $z,\zeta$ [dB]')
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)
legend('$\Sigma$','$\Sigma_{r}$','$\Sigma_{r}^\circ$','$\Omega_0$',...
    'location','SW')
subplot(ax{2,2})
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)

save_fig = 0;
if save_fig
    print('../../Paper/Figures/Ex_Bode','-depsc')
end

%% Plot time domain signals + error bound
% Perfrom simulation using MTF
[rho_NL,zeta_NL]    = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);

% Compute error bound
for i = 1:2
    for k = 1:2
        Hinfnorm(i,k) = norm(sys(i,k)-sysr(i,k),inf);
    end
end
Jbar        = max(max(Hinfnorm));
gammarhou   = norm(sysr(1,1),inf);
gammarhophi = norm(sysr(1,2),inf);
gammazphi   = norm(sysr(2,2),inf);
gammayphi   = norm(sys(1,2),inf);
gammabar    = (1+gammarhou*gamma/(1-gammarhophi*gamma))*...
                (1+gammazphi*gamma/(1-gammayphi*gamma));
nrmU        = signal_norm(u,1);

maxE        = nrmU*Jbar*gammabar;

figure
subplot(211)
plot(tmtf,z_NL,tmtf,zeta_NL)
xlabel('Time [sec]')
ylabel('$\bar z(t), \bar \zeta(t)$')
legend('$\bar z(t)$','$\bar \zeta(t)$')

subplot(212)
plot(tmtf,z_NL-zeta_NL)
xlabel('Time [sec]')
ylabel('$\bar z - \bar \zeta$(t)')

display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL-zeta_NL,1)) ' <= ' num2str(maxE)])

%% Validation test 1 
ftest1 = 5;
nH     = 10;
A      = 1:nH;
A      = A*0+1;

% Define time vector
Tmtf_val1 = 1/ftest1;           % End time
ts        = 1e-4;               % Sampling time
tmtf_val1 = 0:ts:Tmtf_val1-ts;  % Time vector

% Define input vector
u_val1 = tmtf_val1*0;
for i = 1:nH
    u_val1 = u_val1 + A(i)*sin(2*pi*ftest1*i*tmtf_val1);
end
plot_u_val = 0;
if plot_u_val
    plot(tmtf_val1,u_val1)
    xlabel('Time [sec]')
    ylabel('$u(t)$')
end

nrmU_val1       = signal_norm(u_val1,1);
maxE_val1_nrm   = nrmU_val1*Jbar*gammabar;

[~,z_NL_val1]       = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
[~,zeta_NL_val1]    = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
[~,zeta_NL_val01]   = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,@(y)(y*0));

%% Validation test 2
ftest2 = 100;
nH     = 10;
A      = 1:nH;
A      = A*0+1;

% Define time vector
Tmtf_val2 = 1/ftest2;           % End time
ts        = 1e-4;               % Sampling time
tmtf_val2 = 0:ts:Tmtf_val2-ts;  % Time vector

% Define input vector

u_val2 = tmtf_val2*0;
for i = 1:nH
    u_val2 = u_val2 + A(i)*sin(2*pi*ftest2*i*tmtf_val2);
end
plot_u_val = 0;
if plot_u_val
    plot(tmtf,u_val2)
    xlabel('Time [sec]')
    ylabel('$u(t)$')
end

nrmU_val2     = signal_norm(u_val2,1);
maxE_val2_nrm = nrmU_val2*Jbar*gammabar;

[~,z_NL_val2]       = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
[~,zeta_NL_val2]    = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
[~,zeta_NL_val02]   = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,@(y)(y*0));

%% Perform MOR for generic NL models
% Polynomial fit for h \circ \pi
N       = size(omega,1);
phi     = [ones(1,N); omega'; (omega').^2; omega(:,1)'.*omega(:,2)'; ...
    omega(:,1)'.*omega(:,3)';omega(:,1)'.*omega(:,4)';omega(:,1)'.*omega(:,5)';
    omega(:,1)'.*omega(:,6)';omega(:,1)'.*omega(:,7)';
    omega(:,2)'.*omega(:,3)';omega(:,2)'.*omega(:,4)';omega(:,2)'.*omega(:,5)';
    omega(:,2)'.*omega(:,6)';omega(:,2)'.*omega(:,7)';
    omega(:,3)'.*omega(:,4)';omega(:,3)'.*omega(:,5)';omega(:,3)'.*omega(:,6)';
    omega(:,3)'.*omega(:,7)';
    omega(:,4)'.*omega(:,5)';omega(:,4)'.*omega(:,6)';omega(:,4)'.*omega(:,7)';
    omega(:,5)'.*omega(:,6)';omega(:,5)'.*omega(:,7)';
    omega(:,6)'.*omega(:,7)';];
alpha = z_NL'/phi;

% Display a measure for the mismatch between h \circ \pi and the fitted
% function
norm(z_NL-(alpha*phi)')

%% Define reduced-order Wiener model
sysr_gen = ss(sysr.A(1:v{1,1},1:v{1,1}),sysr.B(1:v{1,1},1),sysr.C(1,1:v{1,1}),0);
varphi   = @(omega)(alpha*[ones(1,size(omega,1)); omega'; (omega').^2; omega(:,1)'.*omega(:,2)'; ...
    omega(:,1)'.*omega(:,3)';omega(:,1)'.*omega(:,4)';omega(:,1)'.*omega(:,5)';
    omega(:,1)'.*omega(:,6)';omega(:,1)'.*omega(:,7)';
    omega(:,2)'.*omega(:,3)';omega(:,2)'.*omega(:,4)';omega(:,2)'.*omega(:,5)';
    omega(:,2)'.*omega(:,6)';omega(:,2)'.*omega(:,7)';
    omega(:,3)'.*omega(:,4)';omega(:,3)'.*omega(:,5)';omega(:,3)'.*omega(:,6)';
    omega(:,3)'.*omega(:,7)';
    omega(:,4)'.*omega(:,5)';omega(:,4)'.*omega(:,6)';omega(:,4)'.*omega(:,7)';
    omega(:,5)'.*omega(:,6)';omega(:,5)'.*omega(:,7)';
    omega(:,6)'.*omega(:,7)';]);

%% Simulate Wiener model
zeta_NL_pol      = sim_Wiener(u,tmtf,sysr_gen,varphi)';
zeta_NL_pol_val1 = sim_Wiener(u_val1,tmtf_val1,sysr_gen,varphi)';
zeta_NL_pol_val2 = sim_Wiener(u_val2,tmtf_val2,sysr_gen,varphi)';

%% Plot training results
figure
subplot(211)
plot(tmtf,z_NL,tmtf,zeta_NL_pol,'--')
hold all
plot(tmtf,zeta_NL,':')
xlabel('Time [sec]')
legend('$\bar z$','$\bar \zeta_{Wiener}$','$\bar \zeta_{Lur''e}$')
ylabel('Deflection [m]')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

subplot(212)
plot(tmtf,z_NL-zeta_NL_pol,tmtf,z_NL-zeta_NL,'--')
xlabel('Time [sec]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$')
ylabel('Error [m]')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper/Figures/Ex_Time_Training','-depsc')
end

%% Plot validation results for val1
figure
subplot(211)
plot(tmtf_val1,z_NL_val1,tmtf_val1,zeta_NL_pol_val1,'--')
hold all
plot(tmtf_val1,zeta_NL_val1,':')
xlabel('Time [sec]')
ylabel('Deflection [m]')
legend('$\bar z$','$\bar \zeta_{Wiener}$','$\bar \zeta_{Lur''e}$')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

subplot(212)
plot(tmtf_val1,z_NL_val1-zeta_NL_pol_val1)
hold all
plot(tmtf_val1,z_NL_val1-zeta_NL_val1,'--')
xlabel('Time [sec]')
ylabel('Error [m]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper/Figures/Ex_Time_Val1','-depsc')
end

%% Plot validation results for val2
figure
subplot(211)
plot(tmtf_val2,z_NL_val2,tmtf_val2,zeta_NL_pol_val2,'--')
hold all
plot(tmtf_val2,zeta_NL_val2,':')
xlabel('Time [sec]')
ylabel('Deflection [m]')
legend('$\bar z$','$\bar \zeta_{Wiener}$','$\bar \zeta_{Lur''e}$')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

subplot(212)
plot(tmtf_val2,z_NL_val2-zeta_NL_pol_val2)
hold all
plot(tmtf_val2,z_NL_val2-zeta_NL_val2,'--')
xlabel('Time [sec]')
ylabel('Error [m]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper/Figures/Ex_Time_Val2','-depsc')
end

%% Values for table
display('Training results')
display(['|z|_2 = ' num2str(signal_norm(z_NL,1))])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL-zeta_NL,1)) ' <= ' num2str(maxE)])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL-zeta_NL_pol,1))])
display('Validation 1 results')
display(['|z|_2 = ' num2str(signal_norm(z_NL_val1,1))])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1,1)) ' <= ' num2str(maxE_val1_nrm)])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_pol_val1,1))])
display('Validation 2 results')
display(['|z|_2 = ' num2str(signal_norm(z_NL_val2,1))])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2,1)) ' <= ' num2str(maxE_val2_nrm)])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_pol_val2,1))])