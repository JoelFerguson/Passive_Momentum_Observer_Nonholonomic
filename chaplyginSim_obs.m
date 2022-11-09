clear
clc
% Define macro for extracting matrix entries
matrixIdx = @(A,i,j) A(i,j);

%% Define system model
% Parameters
sys.m = 1;
sys.J = 1;
sys.l = 0.5;
sys.d1 = 0;
sys.d2 = 0;

% System matrices and Hamiltonian
sys.Mr = @(q) [sys.m 0; 0 sys.J+sys.m*sys.l^2];
sys.Gr = @(q) eye(2);
sys.Dr = @(q) diag([sys.d1; sys.d2]);
sys.Sr = @(q,p) [0 p(2)*sys.m*sys.l/(sys.J + sys.m*sys.l^2);
                -p(2)*sys.m*sys.l/(sys.J + sys.m*sys.l^2) 0];
sys.Q = @(q) [1 q(3); 0 1; 0 -q(1)];
sys.Q1 = @(q) matrixIdx(sys.Q(q),1:2,1:2);
sys.Q2 = @(q) matrixIdx(sys.Q(q),3,1:2);
sys.V = @(q) 0;
sys.T = @(q,p) 0.5*p.'*(sys.Mr(q)\p);
sys.H = @(q,p) sys.T(q,p) + sys.V(q);

% Energy gradients
syms q1 q2 q3
syms p1 p2
q_sym = [q1; q2; q3];
p_sym = [p1; p2];
% Kinetic energy
sys.dTdq_int = matlabFunction(simplify(jacobian(sys.T(q_sym,p_sym),q_sym)).','vars',[q_sym; p_sym]);
sys.dTdq = @(q,p) sys.dTdq_int(q(1),q(2),q(3),p(1),p(2));
sys.dHdp = @(q,p) sys.Mr(q)\p;
% Potential energy
sys.dVdq_int = matlabFunction(jacobian(sys.V(q_sym),q_sym).','vars',q_sym);
sys.dVdq = @(q) sys.dVdq_int(q(1),q(2),q(3));
% Total energy
sys.dHdq = @(q,p) sys.dTdq(q,p) + sys.dVdq(q);

% Define open-loop dynamics --- x = [q; p]
sys.dq = @(q,p) sys.Q(q)*sys.dHdp(q,p);
sys.dp = @(q,p,u) -sys.Q(q).'*sys.dHdq(q,p) + (sys.Sr(q,p) - sys.Dr(q))*sys.dHdp(q,p) + sys.Gr(q)*u;
sys.dx = @(x,u) [sys.dq(x(1:3),x(4:5)); sys.dp(x(1:3),x(4:5),u)];

%% Define momentum observer
% Set tuning parameter to adjust rate of convergence
obs.kappa = 2;

% Define intermediate inverse mass to normalise kinetic energy
obs.Mbi = @(q) sys.Q1(q)*(sys.Mr(q)\sys.Q1(q).');
% Compute positive square root of Mbi
obs.tau = @(q) matrixIdx(obs.Mbi(q),1,1) + matrixIdx(obs.Mbi(q),2,2);
obs.s = @(q) sqrt(det(obs.Mbi(q)));
obs.t = @(q) sqrt(obs.tau(q) + 2*obs.s(q));
obs.T = @(q) (obs.Mbi(q)+obs.s(q)*eye(2))/obs.t(q);

% Compute matrices for normalised system
obs.A = @(q) obs.T(q)/(sys.Q1(q).');
obs.B = @(q) sys.Q2(q)*obs.A(q);
obs.D = @(q) obs.A(q)*sys.Dr(q)*obs.A(q).';
obs.G = @(q) obs.A(q)*sys.Gr(q);
% Define transformed momentum vector
obs.p = @(q,pr) obs.A(q)*pr;

% Compute coriollis matrix
obs.S_sym = matlabFunction(sys.Q(q_sym).'*jacobian(obs.A(q_sym)\p_sym,q_sym).' - jacobian(obs.A(q_sym)\p_sym,q_sym)*sys.Q(q_sym),'vars',[q_sym; p_sym]);
obs.S = @(q,p) obs.A(q)*sys.Sr(q,obs.A(q)\p)*obs.A(q).' + obs.A(q)*obs.S_sym(q(1),q(2),q(3),p(1),p(2))*obs.A(q).';
% Define additional symbolic terms for computing Sbar
syms pb1 pb2
p_sym2 = [pb1 pb2].';
% Compute Sbar
obs.SbSYM = matlabFunction(jacobian(obs.S(q_sym,p_sym)*p_sym2,p_sym),'vars',[q_sym; p_sym2]);
obs.Sb = @(q,ph) obs.SbSYM(q(1),q(2),q(3),ph(1),ph(2));

% Define energy function of observer
obs.ph = @(q,xp,phi) xp + phi*q(1:2);
obs.Ho = @(p,ph) 0.5*(p-ph).'*(p-ph);

% Define condition for switching
obs.switchCond = @(q,ph,phi) min(eig(phi*eye(2)*obs.T(q) - 0.5*(obs.Sb(q,ph) + obs.Sb(q,ph).'))) - obs.kappa;
obs.stopEventWrapper = @(t,x) stopEvent(t,x,obs.switchCond,obs.ph);

% Construct observer
obs.dxp = @(q,ph,u,uo,phi) (obs.S(q,ph) - obs.D(q) - phi*eye(2)*obs.T(q))*ph - [obs.T(q) obs.B(q).']*sys.dVdq(q) + obs.G(q)*(u + uo);
obs.dphi = 0;

%% Simulate system
% Define system input
u = @(t,q) [1*sin(10*t);2*cos(4*t)];
% Define observer input
uo = @(q) [0;0];

% Simulation time
sim.t = [0 4];

% Initial conditions for q, p, phat. Note that phat is in transformed
% coordinates
sim.q0 = [0; 0; 0];
sim.p0 = [2; -1];
sim.ph0 = [0 0].';
% compute initial phi to ensure system starts in flow dynamics
sim.phi0 = 0;
while obs.switchCond(sim.q0,sim.ph0,sim.phi0) <= 0
    sim.phi0 = sim.phi0 + obs.kappa;
end
% Construct initial value for observer state xp
sim.xp0 = sim.ph0 - sim.phi0*sim.q0(1:2);

% Construct simulation initial conditions
sim.x0 = [sim.q0; sim.p0; sim.xp0; sim.phi0];

% Construct overall ode --- x = (q,p,xp,phi)
ode = @(t,x) [sys.dx(x(1:5),u(t,x(1:3)));
                obs.dxp(x(1:3),obs.ph(x(1:3),x(6:7),x(8)),u(t,x(1:3)),uo(x(1:3)),x(8));
                obs.dphi];
% stop solver when gain update is required
options = odeset('RelTol',1e-6,'Events',obs.stopEventWrapper);
% Start solver
[t,x] = ode45(ode,sim.t,sim.x0,options);
% Store outputs
res.x = x;
res.t = t;

% Iteratively continue solve until final time point reached
while(t(end) < sim.t(2))
    % Extract initial conditions for next flow phase
    x_end = x(end,:).';
    q_end = x_end(1:3);
    p_end = x_end(4:5);
    xp_end = x_end(6:7);
    phi_end = x_end(8);
    ph_end = obs.ph(q_end,xp_end,phi_end);
    
    % update observer values
    phi_minus = phi_end;
    phi_plus = phi_minus + obs.kappa;
    xp_plus = xp_end + (phi_minus - phi_plus)*q_end(1:2);
    
    % Construct new initial state vector
    sim.x0 = [q_end; p_end; xp_plus; phi_plus];
    
    % Start solver
    [t,x] = ode45(ode,[t(end) sim.t(2)],sim.x0,options);
    
    % Store outputs
    res.x = [res.x; x];
    res.t = [res.t; t];
end

% Unpack results
res.z = res.x(:,1:3);
res.pr = res.x(:,4:5);
res.xp = res.x(:,6:7);
res.phi = res.x(:,8);
% Map configuration vector back to original coordinates
res.q = zeros(size(res.z));
res.H = zeros(length(res.t),1);
res.Ho = zeros(length(res.t),1);
res.p = zeros(size(res.pr));
res.ph = zeros(size(res.pr));
for i=1:length(res.t)
    res.q(i,3) = res.z(i,2);
    res.q(i,1:2) = ([cos(res.z(i,2)) sin(res.z(i,2)); -sin(res.z(i,2)) cos(res.z(i,2))].'*[res.z(i,1); res.z(i,3)]).';
    res.H(i) = sys.H(res.z(i,:).',res.pr(i,:).');
    res.p(i,:) = obs.p(res.z(i,:).',res.pr(i,:).').';
    res.ph(i,:) = obs.ph(res.z(i,:).',res.xp(i,:).',res.phi(i)).';
    res.Ho(i) = obs.Ho(res.p(i,:).',res.ph(i,:).');
end

% figure(1)
% plot(res.q(:,1),res.q(:,2))
% grid on
% axis equal

fig2 = figure(2)
plot(res.t,res.p)
grid on
hold on
plot(res.t,res.ph,'--')
legend('$$p_1$$','$$p_2$$','$$\hat{p}_1$$','$$\hat{p}_2$$','Interpreter','Latex')
xlabel('time (s)')
ylabel('Momentum')

set(findall(fig2,'type','text'),'FontSize',11)
set(findall(fig2,'type','axes'),'FontSize',11)
set(findall(fig2,'type','line'),'linewidth',2.0)

%% Functions
function [value, isTerminal, direction] = stopEvent(t,x,stop_func,ph_func)
    % Stop integration when zero detected    
    isTerminal = 1;
    
    q = x(1:3);
    xp = x(6:7);
    phi = x(8);
    ph = ph_func(q,xp,phi);
    value = stop_func(q,ph,phi);
    
    % detect all zeros
    direction = 0;
end