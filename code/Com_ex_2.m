%% Control 2022
% Computer exercise 2
% Sujet Phodapol 
% Matthew Lock

%% clear
clc;
clear all;
close all;

%% ========== 3.1 Poles, zeros and RGA ============
%% ===== minimum phase =====
%% 3.1.1
% create system with pre-defined function
sysmp = minphase;
[A,B,C,D] = ssdata(sysmp);

% add label for clarity
states = {'tank1' 'tank2' 'tank3' 'tank4'};
inputs = {'Pump1 (u1)' 'Pump2 (u2)'};
outputs = {'Level tank 1 (y1)' 'Level tank 2 (y2)'};

% create state-space model
sys_mimo = ss(A,B,C,D,'statename',states,'inputname',inputs,...
'outputname',outputs);

% find transfer function
G_min = tf(sys_mimo);

% calculate poles & zeros
p11 = pole(G_min(1,1));
z11 = zero(G_min(1,1));

%% 3.1.2 compute poles & zeros of MIMO system
fprintf('=== Poles === \n');
pole(G_min)
fprintf('=== Zeros === \n');
tzero(G_min)

%% 3.1.3 singular value
figure,sigma(G_min)

%% 3.1.4 RGA
rga = G_min.*(inv(G_min))';
rga_0 = evalfr(rga,0);
rga_wc = evalfr(rga,1i*0.1);

%% 3.1.5 step
figure, step(G_min)

%% ===== non-minimum phase =====
%% 3.1.1
% create system with pre-defined function
sysnmp = nonminphase;
[A,B,C,D] = ssdata(sysnmp);

% add label for clarity
states = {'tank1' 'tank2' 'tank3' 'tank4'};
inputs = {'Pump1 (u1)' 'Pump2 (u2)'};
outputs = {'Level tank 1 (y1)' 'Level tank 2 (y2)'};

% create state-space model
sys_mimo = ss(A,B,C,D,'statename',states,'inputname',inputs,...
'outputname',outputs);

% find transfer function
G_nmin = tf(sys_mimo);

% calculate poles & zeros
p11 = pole(G_nmin(1,1));
z11 = zero(G_nmin(1,1));

%% 3.1.2 compute poles & zeros of MIMO system
fprintf('=== Poles === \n');
pole(G_nmin)
fprintf('=== Zeros === \n');
tzero(G_nmin)

%% 3.1.3 singular value
figure,sigma(G_nmin)

%% 3.1.4 RGA
rga = G_nmin.*(inv(G_nmin))';
rga_0 = evalfr(rga,0);
rga_wc = evalfr(rga,1i*0.02);

%% 3.1.5 step
figure,step(G_nmin)
%% ========== 3.2 Decentralized Control ============
%% ===== minimum phase =====
%% 3.2.1 Decentralized controller
%% f1
phi_m = pi/3;
wc = 0.1;
T1 = tan(-pi/2 + phi_m - angle(evalfr(G_min(1,1),1i*wc)))/wc;
f1 = 1+1/(s*T1);
[m,phase] = bode(f1*G_min(1,1),wc);
K1 = 1/m;
f1 = K1*f1;
figure,margin(f1*G_min(1,1))

%% f2
T2 = tan(-pi/2 + phi_m - angle(evalfr(G_min(2,2),1i*wc)))/wc;
f2 = 1+1/(s*T2);
[m,phase] = bode(f2*G_min(2,2),wc);
K2 = 1/m;
f2 = K2*f2;
figure,margin(f2*G_min(2,2))

% create controller matrix
F = [f1 0;0 f2];
G = G_min;  

%% 3.2.2 Decentralized controller
S = inv(eye(2) + G_min * F);
T = inv(eye(2) + G_min * F)*G_min*F;

figure;
sigma(S);
hold;
sigma(T);

L = G*F;

figure;
set(gcf, 'Position', get(0, 'Screensize'));

h=subplot(2,1,1);
margin(L(1,1));
grid;
text('FontSize', 24);
myAxes=findobj(h,'Type','Axes');
exportgraphics(myAxes,'margin_L11.pdf');

h=subplot(2,1,2);
margin(L(2,2));
grid;
myAxes=findobj(h,'Type','Axes');
exportgraphics(myAxes,'margin_L22.pdf');

%% ===== non-minimum phase =====
%% 3.2.1 Decentralized controller
%% f1
phi_m = pi/3;
wc = 0.02;
T1 = tan(-pi/2 + phi_m - angle(evalfr(G_nmin(1,2),1i*wc)))/wc;
f1 = 1+1/(s*T1);
[m,phase] = bode(f1*G_nmin(1,2),wc);
K1 = 1/m;
f1 = K1*f1;
figure,margin(f1*G_nmin(1,2))

%% f2
T2 = tan(-pi/2 + phi_m - angle(evalfr(G_nmin(2,1),1i*wc)))/wc;
f2 = 1+1/(s*T2);
[m,phase] = bode(f2*G_nmin(2,1),wc);
K2 = 1/m;
f2 = K2*f2;
figure,margin(f2*G_nmin(2,1))

% create controller matrix
F = [0 f2;f1 0];
G = G_nmin;

%% 3.2.2 Decentralized controller

L = minreal(G*F);
S = inv(eye(2) + L);
T = (eye(2) + L)^-1*L;
T = minreal(T);

figure;
sigma(S);
hold;
sigma(T);

% figure;
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% h=subplot(2,1,1);
% margin(L(1,1));
% grid;
% text('FontSize', 24);
% myAxes=findobj(h,'Type','Axes');
% exportgraphics(myAxes,'margin_L11.pdf');
% 
% h=subplot(2,1,2);
% margin(L(2,2));
% grid;
% myAxes=findobj(h,'Type','Axes');
% exportgraphics(myAxes,'margin_L22.pdf');

%%
figure('Name', 'Bode');
margin(L(1, 1));
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf, 'Position', get(0, 'Screensize'));
grid;
saveas(gcf, 'L11_margin.png');

figure('Name', 'Bode')
margin(L(2, 2));
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf, 'Position', get(0, 'Screensize'));
grid;
saveas(gcf, 'L22_margin.png')

%% 3.2.3 Simulink
closedloop

%% plot simulink

y1 = yout.Data(:,1);
y2 = yout.Data(:,2);

u1 = uout.Data(:,1);
u2 = uout.Data(:,2);

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(tout,y1)
hold on;
plot(tout,y2)
legend('y1','y2');
title('System output step response');
xlabel('Time (s)');
ylabel('Output');
grid;
set(findall(gcf,'-property','FontSize'),'FontSize',24)
% saveas(gcf, 'nonmin_phase_y.pdf')

figure;
set(gcf, 'Position', get(0, 'Screensize'));
plot(tout,u1)
hold on;
plot(tout,u2)
legend('u1','u2');
title('Control output step response');
xlabel('Time (s)');
ylabel('Output');
grid;
set(findall(gcf,'-property','FontSize'),'FontSize',24)
% saveas(gcf, 'nonmin_phase_u.pdf')




%% ========= EXCESSS ==========
%% test
Gs = 1/(s+1)*[1 -1;1.1 -1];
rgas = Gs.*(1/Gs)';
wc = 0;
val = evalfr(rgas,j*wc)

%% test 2
G_min = tf([1 2],[3 4 5]) %(as an example))
w = 0; % for example 3 radians/sec
val = evalfr(G_min,j*w)

%% Plot
s = tf('s');
sys = 1/(s+5);
figure;
bode(sys)
set(findall(gcf,'type','line'),'linewidth',2);
grid on