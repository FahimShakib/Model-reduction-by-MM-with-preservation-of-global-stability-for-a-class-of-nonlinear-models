%% Description
% This files contains the example in the paper [] on model order
% reduction by moment matching for convergent Lur'e-type models.
% Version 2 includes a comparison with balanced truncation.
% Version 3 includes tangential moment matching

% Author: Fahim Shakib
% Date:   Feb. 07, 2023
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

Ne  = 200;          % Number of beam elements

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

% Use double variables for tangential moment matching (Refered to as MIMO MM)
% Only use interpolation points Omega0{1,1}
SS = S{1,1};
LL = [L{1,1}; L{2,1}];

%% Input to MTF 
% MTF algorithm is used to simulate the Lur'e-type system below
% See [Pavlov A, Hunnekens BG, Wouw N, Nijmeijer H. Steady-state performance
% optimization for nonlinear control systems of Lur’e type. Automatica. 
% 2013 Jul 1;49(7):2087-97.] for the MTF algorithm

% MTF Settings
Tmtf      = 5;              % End time
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
load_hinf = 1;
if load_hinf
    load('Data/20230221 - Results paper','nrm_sys')
else
    nrm_sys = hinfnorm(sys(1,2));
end
display(['Max. \gamma^\star = ' num2str(1/nrm_sys)])

gamma   = 8.6e4;                                % Spring stiffness
NLfnc   = @(y)onesidedspring(y,gamma);          % Define one-sided spring

display(['Verify the convergence property:'])
display(['gamma*|Sys|_infty = ' num2str(gamma*nrm_sys) ' < 1'])

%% Generate data NL System
load_y_NL = 1;
% Perfrom simulation using MTF
if load_y_NL
    load('Data/20230221 - Results paper','y_NL','z_NL','y_NL0','z_NL0','t_z_NL')
else
    % Simulate model
    tic
    [y_NL,z_NL]    = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);
    t_z_NL = toc;
    % Simulation without nonlinearity (the nonlinearity always outputs 0),
    % just to see the effect of the nonlinearity on the steady-state
    % response.
    [y_NL0,z_NL0]  = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,@(y)(y*0));
end
% Plot both simulation responses
plot_z_NL = 1;
if plot_z_NL
    figure
    plot(tmtf,z_NL,tmtf,z_NL0)
    xlabel('Time [sec]')
    ylabel('$\bar z$(t)')
    legend('With nonlinearity','Without nonlinearity')
end

%% SISO MM: Perform MOR - Solve Sylvester equations to compute CPi
% by Lyap (AX + XB + C = 0)
for i = 1:2
    for k = 1:2
        Pi{i,k}  = lyap(sys.A,-S{i,k},sys.B(:,k)*L{i,k});
        CPi{i,k} = sys.C(i,:)*Pi{i,k};
        H{i,k}   = CPi{i,k};
    end
end

%% MIMO MM: Perform tangential MOR - Solve Sylvester equations to compute CPi
PiMIMO    = lyap(sys.A,-SS,sys.B*LL);
CPiMIMO   = sys.C*PiMIMO;
HMIMO     = CPiMIMO;

%% MIMO MM: Initialize to dominant poles of A
sigmaA  = eig(sys.A);
pp      = [sigmaA(end-v{1,1}+1:end);];

GinitMIMO       = place(SS',LL',pp)';
FinitMIMO       = SS-GinitMIMO*LL;
sysrinitMIMO    = ss(FinitMIMO,GinitMIMO,HMIMO,0);

display(['Verify the convergence property: gamma*|Sys_red|_infty = '])
display([num2str(gamma*norm(sysrinitMIMO(1,2),inf)) ' < 1'])

bodeplot = 0;
if bodeplot
    figure
    bode(sys)
    hold all
    bode(sysrinitMIMO)
    legend('Full-order model','Reduced (initial) model by MIMO MM')
end

% Initial model found by MIMO MM is not convergent, which is expected as
% there are no constraints that ensure that this model is convergent.

%% SISO MM: Find initial model - The approach proposed in our paper
% Define sdp variables
clear P
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

display(['Verify the convergence property:'])
display(['gamma*|Sys_red|_infty = ' num2str(gamma*norm(sysrinit(1,2),inf)) ' < 1'])

%% Inspect FRF for match
% Display the FRF at interpolation frequencies
for i = 1:2
    for k = 1:2
        disp(['Sys0(' num2str(i) ',' num2str(k) ')-Sys_red(' num2str(i) ',' num2str(k) '): '...
            num2str(abs(freqresp(sys(i,k)-sysrinit(i,k),2*pi*Omega0{i,k})))])
    end
end

%% Solve optimization problem using a penalty function approach
% Select frequencies to minimize cost function on
wlist = unique([logspace(-2,5,1000) Omega0{1,1} Omega0{1,2} Omega0{2,1} Omega0{2,2}])*2*pi;

% Define cost function
[mag,~,wout] = bode(sys,wlist);
options      = sdpsettings('solver','mosek','verbose',0);
J_SISO       = @(theta,imp)cost_fnc_FRF_penalty(theta,H,S,L,mag,wout,v,gamma,imp);
J_SISO_unc   = @(theta)cost_fnc_FRF(theta,H,S,L,mag,wout,v,gamma);
J_MIMO       = @(theta)cost_fnc_FRF_penalty_MIMO(theta,HMIMO,SS,LL,mag,wout,v,gamma);

clear options
% Optimization settings
options                         = optimoptions('fminunc','Display','iter');
options.MaxFunctionEvaluations  = 10000;
options.OptimalityTolerance     = 1e-16;
% Perform optimization

% Approached proposed in the paper. The constrained optization problem is
% solved using a penalty-function approach.
tic
thetaopt_SISO                   = fminunc(J_SISO,thetainit,options,0);
t_hinf = toc;
% Same approached as proposed in the paper, but with a different
% implementataion of the constraints. With this different implementation,
% it takes longer to solve the problem.
tic
% thetaopt_SISO                   = fminunc(J_SISO,thetainit,options,1);
t_lmi = toc;
% Same approached as proposed in the paper, but without constraints.
% tic
% thetaopt_SISO_unc               = fminunc(J_SISO_unc,thetainit,options);
% t_unc = toc;
% MIMO moment matching approach using tangential moment matching and a
% penalty function to enforce constraints
tic
thetaopt_MIMO                   = fminunc(J_MIMO,GinitMIMO,options);
t_MIMO = toc;

% Construct reduced-order LTI model
[~,sysr]     = J_SISO(thetaopt_SISO,0);
% [~,sysr_unc] = J_SISO_unc(thetaopt_SISO_unc);
[~,sysrMIMO] = J_MIMO(thetaopt_MIMO);

display(['Verify the convergence property:'])
display(['gamma*|Sys_red|_infty = ' num2str(gamma*norm(sysr(1,2),inf)) ' < 1'])
display(['gamma*|SysMIMO_red|_infty = ' num2str(gamma*norm(sysrMIMO(1,2),inf)) ' < 1'])

display(['Spent time in optimizer:'])
display(['Sys_red = ' num2str(t_hinf)])
display(['SysMIMO_red = ' num2str(t_MIMO)])

%% Perform balanced truncation
% For comparison: order for reduction with balanced truncation
n_bt = 7;
n_bt2 = 4*n_bt;

hsing       = hankelsv(sys);
hsing_err   = 2*sum(hsing(n_bt+1:end));
hsing_err2  = 2*sum(hsing(n_bt2+1:end));
sys_bt      = balred(sys,n_bt);
sys_bt2     = balred(sys,n_bt2);

%% Convergence preserved according to a priori error bounds?
ii = 1;
conv_bt = 0;
while ~conv_bt
    if (nrm_sys+2*sum(hsing(ii+1:end)))*gamma < 1
        conv_bt = 1;
    end
    ii = ii + 1;
end
n_bt_min_a_priori = ii;
display(['Minimal order BT (a priori) = ' num2str(n_bt_min_a_priori)])

% Convergence preserved (a posteriori check)
display(['Verify the convergence property a posteriori:'])
display(['gamma*|Sys_bt|_infty = ' num2str(gamma*norm(sys_bt(1,2),inf)) ' < 1'])
display(['gamma*|Sys_bt2|_infty = ' num2str(gamma*norm(sys_bt2(1,2),inf)) ' < 1'])

%% Check for which reduction orders the reduced-order model is convergent 
% for BT using an a posteriori check
check_convergence_BT = 0;
if check_convergence_BT
    for ii = 0:50
        sys_bt_tmp = balred(sys,ii);
        conv_bt_check(ii+1) = gamma*norm(sys_bt_tmp(1,2),inf) < 1;
    end
    figure
    plot(0:50,conv_bt_check,'x')
    ylim([0 1.1])
    ylabel('Convergent')
    xlabel('Order of reduction model with balanced truncation')
end

%% Plot final bode (all ROMs)
% Compute FRF of full-order LTI model
mag = bode(sys,wlist);
Mag = squeeze(20*log10(mag));

% Compute FRF of final reduced-order LTI model (SISO MM)
magr = bode(sysr,wlist);
Magr = squeeze(20*log10(magr));

% Compute FRF of final reduced-order LTI model (MIMO MM)
magrMIMO = bode(sysrMIMO,wlist);
MagrMIMO = squeeze(20*log10(magrMIMO));

% Compute FRF of BT7 LTI model
mag_bt = bode(sys_bt,wlist);
Mag_bt = squeeze(20*log10(mag_bt));

% Compute FRF of BT28 LTI model
mag_bt2 = bode(sys_bt2,wlist);
Mag_bt2 = squeeze(20*log10(mag_bt2));

fig = figure;
fig.Position = [680 678 1172 420];
kk = 1;
for i = 1:2
    for k = 1:2
        ax{i,k} = subplot(2,2,kk);
        kk = kk + 1;
        semilogx(wlist/2/pi,squeeze(Mag(i,k,:)),'LineWidth',2);hold all
        semilogx(wlist/2/pi,squeeze(Magr(i,k,:)),'LineWidth',2,'linestyle','--');
        semilogx(wlist/2/pi,squeeze(MagrMIMO(i,k,:)),'LineWidth',2,'linestyle',':');
        semilogx(wlist/2/pi,squeeze(Mag_bt(i,k,:)),'LineWidth',2,'linestyle','-');
        semilogx(wlist/2/pi,squeeze(Mag_bt2(i,k,:)),'LineWidth',2,'linestyle','-.');
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
legend('$\Sigma$','$\Sigma_{SISO-MM}$','$\Sigma_{MIMO MM}$',...
    '$\Sigma_{BT7}$','$\Sigma_{BT28}$','$\Omega_0$','location','SW')
subplot(ax{2,2})
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Bode_RESP','-depsc')
end

%% Plot tangential directions
% Compute FRF of full-order LTI model
systang = ss(sys.a,[sys.b sys.b*LL(:,1)],sys.c,0);
sysrinitMIMOtang = ss(sysrinitMIMO.a,[sysrinitMIMO.b sysrinitMIMO.b*LL(:,1)],sysrinitMIMO.c,0);

magtang = bode(systang,wlist);
Magtang = squeeze(20*log10(magtang));

% Compute FRF of final reduced-order LTI model (MIMO MM)
magrinitMIMOtang = bode(sysrinitMIMOtang,wlist);
MagrinitMIMOtang = squeeze(20*log10(magrinitMIMOtang));

fig = figure;
fig.Position = [680 678 1172 420];
kk = 1;
for i = 1:2
    for k = 1:3
        ax{i,k} = subplot(2,3,kk);
        kk = kk + 1;
        semilogx(wlist/2/pi,squeeze(Magtang(i,k,:)),'LineWidth',2);hold all
        semilogx(wlist/2/pi,squeeze(MagrinitMIMOtang(i,k,:)),'LineWidth',2,'linestyle','--');
        magIP = squeeze(20*log10(bode(systang(i,k),2*pi*Omega0{1,1})));
        semilogx(Omega0{1,1},magIP,'kx','MarkerSize', 12)
    end
end
linkaxes([ax{1},ax{2},ax{3},ax{4},ax{5},ax{6}],'xy')
xlim([min(wlist) max(wlist)]/2/pi)
xlim([1 1e3])
ylim([-150 -70])

subplot(ax{1,1})
ylabel('Output $y,\rho$ [dB]')
title('Input $u$')
set(gca,'fontsize', 12)
subplot(ax{1,2})
title('Input $\varphi(y),\varphi(\rho)$')
set(gca,'fontsize', 12)
subplot(ax{1,3})
title('Input $u + \varphi(y),u + \varphi(\rho)$')
set(gca,'fontsize', 12)
subplot(ax{2,1})
ylabel('Output $z,\zeta$ [dB]')
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)
legend('$\Sigma$','$\Sigma_{MIMO-MM}^\circ$',...
    '$\Omega_0$','location','SW')
subplot(ax{2,2})
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)
subplot(ax{2,3})
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Bode_tang_RESP','-depsc')
end

%% Plot Bode (paper)
% Compute FRF of initial reduced-order LTI model
magrinit = bode(sysrinit,wlist);
Magrinit = squeeze(20*log10(magrinit));

% Compute FRF of full-order LTI model
mag = bode(sys,wlist);
Mag = squeeze(20*log10(mag));

% Compute FRF of final reduced-order LTI model
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
ylim([-250 -70])

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
legend('$\Sigma$','$\Sigma_{r}$','$\Sigma_{r}^\circ$',...
    '$\Omega_0$','location','SW')
subplot(ax{2,2})
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Bode','-depsc')
end

%% Compute time domain signals + error bound
% Perfrom simulation using MTF
tic
[rho_NL,zeta_NL]        = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);
t_zeta_NL = toc;
tic
[rho_NL,zeta_NL_MIMO]   = MTF(sysrMIMO.A,sysrMIMO.B(:,2),sysrMIMO.B(:,1),sysrMIMO.C(1,:),sysrMIMO.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);
t_zeta_NL_MIMO = toc;
tic
[rho_bt,zeta_bt]        = MTF(sys_bt.A,sys_bt.B(:,2),sys_bt.B(:,1),sys_bt.C(1,:),sys_bt.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);
t_zeta_bt = toc;
tic
[rho_bt2,zeta_bt2]      = MTF(sys_bt2.A,sys_bt2.B(:,2),sys_bt2.B(:,1),sys_bt2.C(1,:),sys_bt2.C(2,:),max_iter,tol,1,Tmtf,length(tmtf),u,NLfnc);
t_zeta_bt2 = toc;

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

%% Pot time domain signals + error bound
fig = figure;
fig.Position = [680 678 1172 420];
p1 = plot(tmtf,z_NL-zeta_NL_MIMO); hold all
p2 = plot(tmtf,z_NL-zeta_bt);
p3 = plot(tmtf,z_NL-zeta_NL);
p4 = plot(tmtf,z_NL-zeta_bt2);
xlabel('Time [sec]')
ylabel('Error')
legend([p3 p1 p2 p4],'$\bar z - \bar \zeta_{SISO-MM}$(t)',...
       '$\bar z - \bar \zeta_{MIMO-MM}$(t)',...
       '$\bar z - \bar \zeta_{BT7}$(t)',...
       '$\bar z - \bar \zeta_{BT28}$(t)')
   
set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL-zeta_NL,1))])
display(['|z-zeta_MIMO|_2 = ' num2str(signal_norm(z_NL-zeta_NL_MIMO,1))])
display(['|z-zeta_bt|_2 = ' num2str(signal_norm(z_NL-zeta_bt,1))])
display(['|z-zeta_bt2|_2 = ' num2str(signal_norm(z_NL-zeta_bt2,1))])

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Time_Training_RESP','-depsc')
end

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

tic
[~,z_NL_val1]         = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
t_z_NL_val1 = toc;
tic
[~,zeta_NL_val1]      = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
t_zeta_NL_val1 = toc;
tic
[~,zeta_NL_MIMO_val1] = MTF(sysrMIMO.A,sysrMIMO.B(:,2),sysrMIMO.B(:,1),sysrMIMO.C(1,:),sysrMIMO.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
t_zeta_NL_MIMO_val1 = toc;
tic
[~,zeta_NL_val1_bt]   = MTF(sys_bt.A,sys_bt.B(:,2),sys_bt.B(:,1),sys_bt.C(1,:),sys_bt.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
t_zeta_NL_val1_bt = toc;
tic
[~,zeta_NL_val1_bt2]  = MTF(sys_bt2.A,sys_bt2.B(:,2),sys_bt2.B(:,1),sys_bt2.C(1,:),sys_bt2.C(2,:),max_iter,tol,1,Tmtf_val1,length(tmtf_val1),u_val1,NLfnc);
t_zeta_NL_val1_bt2 = toc;

%% Plot results Validation test 1
fig = figure;
fig.Position = [680 678 1172 420];
p1 = plot(tmtf_val1,z_NL_val1-zeta_NL_MIMO_val1); hold all
p2 = plot(tmtf_val1,z_NL_val1-zeta_NL_val1_bt);
p3 = plot(tmtf_val1,z_NL_val1-zeta_NL_val1);
p4 = plot(tmtf_val1,z_NL_val1-zeta_NL_val1_bt2);
xlabel('Time [sec]')
ylabel('Error')
legend([p3 p1 p2 p4],'$\bar z - \bar \zeta_{SISO-MM}$(t)',...
       '$\bar z - \bar \zeta_{MIMO-MM}$(t)',...
       '$\bar z - \bar \zeta_{BT7}$(t)',...
       '$\bar z - \bar \zeta_{BT28}$(t)',...
       'location','NW')
   
set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

display(['Validation 1: |z-zeta|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1,1))])
display(['Validation 1: |z-zeta_MIMO|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_MIMO_val1,1))])
display(['Validation 1: |z-zeta_bt|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1_bt,1))])
display(['Validation 1: |z-zeta_bt2|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1_bt2,1))])

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Time_Val1_RESP','-depsc')
end

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

tic
[~,z_NL_val2]         = MTF(sys.A,sys.B(:,2),sys.B(:,1),sys.C(1,:),sys.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
t_z_NL_val2 = toc;
tic
[~,zeta_NL_val2]      = MTF(sysr.A,sysr.B(:,2),sysr.B(:,1),sysr.C(1,:),sysr.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
t_zeta_NL_val2 = toc;
tic
[~,zeta_NL_MIMO_val2] = MTF(sysrMIMO.A,sysrMIMO.B(:,2),sysrMIMO.B(:,1),sysrMIMO.C(1,:),sysrMIMO.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
t_zeta_NL_MIMO_val2 = toc;
tic
[~,zeta_NL_val2_bt]   = MTF(sys_bt.A,sys_bt.B(:,2),sys_bt.B(:,1),sys_bt.C(1,:),sys_bt.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
t_zeta_NL_val2_bt = toc;
tic
[~,zeta_NL_val2_bt2]  = MTF(sys_bt2.A,sys_bt2.B(:,2),sys_bt2.B(:,1),sys_bt2.C(1,:),sys_bt2.C(2,:),max_iter,tol,1,Tmtf_val2,length(tmtf_val2),u_val2,NLfnc);
t_zeta_NL_val2_bt2 = toc;

%% Plot results Validation test 2
fig = figure;
fig.Position = [680 678 1172 420];
p1 = plot(tmtf_val2,z_NL_val2-zeta_NL_MIMO_val2); hold all
p2 = plot(tmtf_val2,z_NL_val2-zeta_NL_val2_bt);
p3 = plot(tmtf_val2,z_NL_val2-zeta_NL_val2);
p4 = plot(tmtf_val2,z_NL_val2-zeta_NL_val2_bt2);
xlabel('Time [sec]')
ylabel('Error')
legend([p3 p1 p2 p4],'$\bar z - \bar \zeta_{SISO-MM}$(t)',...
       '$\bar z - \bar \zeta_{MIMO-MM}$(t)',...
       '$\bar z - \bar \zeta_{BT7}$(t)',...
       '$\bar z - \bar \zeta_{BT28}$(t)',...
       'location','NW')
   
set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

display(['Validation 2: |z-zeta|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2,1))])
display(['Validation 2: |z-zeta_MIMO|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_MIMO_val2,1))])
display(['Validation 2: |z-zeta_bt|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2_bt,1))])
display(['Validation 2: |z-zeta_bt2|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2_bt2,1))])

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Time_Val2_RESP','-depsc')
end

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
tic
zeta_NL_pol      = sim_Wiener(u,tmtf,sysr_gen,varphi)';
t_zeta_NL_pol = toc;
tic
zeta_NL_pol_val1 = sim_Wiener(u_val1,tmtf_val1,sysr_gen,varphi)';
t_zeta_NL_pol_val1 = toc;
tic
zeta_NL_pol_val2 = sim_Wiener(u_val2,tmtf_val2,sysr_gen,varphi)';
t_zeta_NL_pol_val2 = toc;

%% Plot training results (all ROMs)
figure
subplot(211)
plot(tmtf,z_NL)
hold all
plot(tmtf,zeta_NL,':')
plot(tmtf,zeta_bt,'-.')
xlabel('Time [sec]')
legend('$\bar z$','$\bar \zeta_{Lur''e}$',...
    '$\bar \zeta_{Lur''e,bt}$')
ylabel('Deflection [m]')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

subplot(212)
plot(tmtf,z_NL-zeta_NL,'--')
hold all
plot(tmtf,z_NL-zeta_bt,':.')
xlabel('Time [sec]')
legend('$\bar z - \bar \zeta_{Lur''e}$',...
    '$\bar z - \bar \zeta_{Lur''e,bt}$')
ylabel('Error [m]')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper/Figures/Ex_Time_Training','-depsc')
end

%% Plot training results (paper)
figure
co = get(gca,'colororder');

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
plot(tmtf,z_NL-zeta_NL_pol,'--','color',co(2,:))
hold all
plot(tmtf,z_NL-zeta_NL,':','color',co(3,:))
xlabel('Time [sec]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$')
ylabel('Error [m]')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Time_Training','-depsc')
end

%% Plot validation results for val1 (paper)
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
plot(tmtf_val1,z_NL_val1-zeta_NL_pol_val1,'--','color',co(2,:))
hold all
plot(tmtf_val1,z_NL_val1-zeta_NL_val1,':','color',co(3,:))
xlabel('Time [sec]')
ylabel('Error [m]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Time_Val1','-depsc')
end

%% Plot validation results for val2 (all ROMs)
figure
subplot(211)
plot(tmtf_val2,z_NL_val2,tmtf_val2,zeta_NL_pol_val2,'--')
hold all
plot(tmtf_val2,zeta_NL_val2,':',tmtf_val2,zeta_NL_val2_bt,'-.')
xlabel('Time [sec]')
ylabel('Deflection [m]')
legend('$\bar z$','$\bar \zeta_{Wiener}$','$\bar \zeta_{Lur''e}$',...
    '$\bar \zeta_{Lur''e,bt}$')

set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

subplot(212)
plot(tmtf_val2,z_NL_val2-zeta_NL_pol_val2)
hold all
plot(tmtf_val2,z_NL_val2-zeta_NL_val2,'--',tmtf_val2,...
    z_NL_val2-zeta_NL_val2_bt,'-.')
xlabel('Time [sec]')
ylabel('Error [m]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$',...
    '$\bar z - \bar \zeta_{Lur''e,bt}$')


set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper/Figures/Ex_Time_Val2','-depsc')
end

%% Plot validation results for val2 (paper)
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
plot(tmtf_val2,z_NL_val2-zeta_NL_pol_val2,'--','color',co(2,:))
hold all
plot(tmtf_val2,z_NL_val2-zeta_NL_val2,':','color',co(3,:))
xlabel('Time [sec]')
ylabel('Error [m]')
legend('$\bar z - \bar \zeta_{Wiener}$','$\bar z - \bar \zeta_{Lur''e}$')


set(gca,'FontSize',14)
set(findall(gcf,'type','line'),'linewidth',2)

save_fig = 0;
if save_fig
    print('../../Paper_V6/Figures/Ex_Time_Val2','-depsc')
end

%% Values for table
display('Training results')
display(['|z|_2 = ' num2str(signal_norm(z_NL,1))])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL-zeta_NL,1)) ' <= ' num2str(maxE)])
display(['|z-zeta-pol|_2 = ' num2str(signal_norm(z_NL-zeta_NL_pol,1))])
display(['|z-zeta-MIMO|_2 = ' num2str(signal_norm(z_NL-zeta_NL_MIMO,1))])
display(['|z-zeta-BT7|_2 = ' num2str(signal_norm(z_NL-zeta_bt,1))])
display(['|z-zeta-BT28|_2 = ' num2str(signal_norm(z_NL-zeta_bt2,1))])

% Computation times
display('Computation times Training')
display(['FOM = ' num2str(t_z_NL) ' sec'])
display(['SISO MM = ' num2str(t_zeta_NL) ' sec'])
display(['Pol = ' num2str(t_zeta_NL_pol) ' sec'])
display(['MIMO MM = ' num2str(t_zeta_NL_MIMO) ' sec'])
display(['BT7 = ' num2str(t_zeta_bt) ' sec'])
display(['BT28 = ' num2str(t_zeta_bt2) ' sec'])


display('Validation 1 results')
display(['|z|_2 = ' num2str(signal_norm(z_NL_val1,1))])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1,1)) ' <= ' num2str(maxE_val1_nrm)])
display(['|z-zeta-pol|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_pol_val1,1))])
display(['|z-zeta-MIMO|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_MIMO_val1,1))])
display(['|z-zeta-BT7|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1_bt,1))])
display(['|z-zeta-BT28|_2 = ' num2str(signal_norm(z_NL_val1-zeta_NL_val1_bt2,1))])

% Computation times
display('Computation times Validation 1')
display(['FOM = ' num2str(t_z_NL_val1) ' sec'])
display(['SISO MM = ' num2str(t_zeta_NL_val1) ' sec'])
display(['Pol = ' num2str(t_zeta_NL_pol_val1) ' sec'])
display(['MIMO MM = ' num2str(t_zeta_NL_MIMO_val1) ' sec'])
display(['BT7 = ' num2str(t_zeta_NL_val1_bt) ' sec'])
display(['BT28 = ' num2str(t_zeta_NL_val1_bt2) ' sec'])

display('Validation 2 results')
display(['|z|_2 = ' num2str(signal_norm(z_NL_val2,1))])
display(['|z-zeta|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2,1)) ' <= ' num2str(maxE_val2_nrm)])
display(['|z-zeta-pol|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_pol_val2,1))])
display(['|z-zeta-MIMO|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_MIMO_val2,1))])
display(['|z-zeta-BT7|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2_bt,1))])
display(['|z-zeta-BT28|_2 = ' num2str(signal_norm(z_NL_val2-zeta_NL_val2_bt2,1))])

% Computation times
display('Computation times Validation 2')
display(['FOM = ' num2str(t_z_NL_val2) ' sec'])
display(['SISO MM = ' num2str(t_zeta_NL_val2) ' sec'])
display(['Pol = ' num2str(t_zeta_NL_pol_val2) ' sec'])
display(['MIMO MM = ' num2str(t_zeta_NL_MIMO_val2) ' sec'])
display(['BT7 = ' num2str(t_zeta_NL_val2_bt) ' sec'])
display(['BT28 = ' num2str(t_zeta_NL_val2_bt2) ' sec'])
