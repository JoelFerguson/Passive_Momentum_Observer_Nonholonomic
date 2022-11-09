clear
clc

load('obs_kappa_05.mat')
res1.t = res.t;
res1.Hp = res.Ho;
res1.phi = res.phi;

load('obs_kappa_1.mat')
res3.t = res.t;
res3.Hp = res.Ho;
res3.phi = res.phi;

load('obs_kappa_2.mat')
res4.t = res.t;
res4.Hp = res.Ho;
res4.phi = res.phi;

% plot observer energy
fig1 = figure(1)
subplot(2,1,1)
plot(res1.t,log(res1.Hp),'-',res3.t,log(res3.Hp),'--',res4.t,log(res4.Hp),':')
grid on
xlabel('time (s)')
ylabel('Observer energy')
legend({'$\ln H_o, \kappa = 0.5$','$\ln H_o, \kappa = 1$','$\ln H_o, \kappa = 2$'},'Interpreter', 'latex')

subplot(2,1,2)
plot([0;res1.t],[0.5;res1.phi],'-',[0;res3.t],[1;res3.phi],'--',[0;res4.t],[2;res4.phi],':')
grid on
xlabel('time (s)')
ylabel('\phi(t)')
legend({'$\phi, \kappa = 0.5$','$\phi, \kappa = 1$','$\phi, \kappa = 2$'},'Interpreter', 'latex')

set(findall(fig1,'type','text'),'FontSize',11)
set(findall(fig1,'type','axes'),'FontSize',11)
set(findall(fig1,'type','line'),'linewidth',2.0)