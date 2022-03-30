%% Numerical analysis Homework 2

%% Problem 1
clc; clear variables; close all;


kgm = 0.026; %Growth rate
pmax = 12000; %Max. sustainable pop.
dt = 10; %time step
t = 1950:dt:2050; % time vector
IV = 2560; %Initial value

%Defining ODE for dpdt:
dpdt = @(p,t) 0.026*(1-p/pmax)*p;
dpdt2 = @(t,p) 0.026*(1-p/pmax)*p;

%  (a) 
%  Forward euler estimation:
yFE = odeFE(dpdt,IV,t);

%  (b)
%  Backward euler estimation
yBE = odeBE(dpdt,IV,t);

%  (c) 
%  Runge Kutta 4th orde\oder estimation
yRK4 = odeRK4(dpdt,IV,t);

%  (d)
%  Matlabs ode45:
[t_ode45,y_ode45] = ode45(dpdt2,t,IV);

%  (e)
%first defining analytical solution:
p0 = 2560;
AS = @(t) p0 * (pmax ./ (p0 + (pmax - p0)*exp(-kgm*(t-1950))));
ASval = AS(t);

%Actual world population:
Year = 1950:5:2000;
Pop = [2560, 2780, 3040, 3350, 3710, 4090, 4450, 4850, 5280, 5690, 6080];

%  plotting results:
plot(t,yFE)
hold on
plot(t,yBE)
plot(t,yRK4)
plot(t_ode45,y_ode45);
plot(t,ASval)
plot(Year, Pop,'r*','MarkerSize',5)
grid
legend('Forward Euler','Backward Euler','Runge Kutta 4th order', 'ode45','Analytical solution', 'Actual Pop','Location','Northwest')
title('Simulation of world population using various methods')
xlabel('Year')
ylabel('population [million people]')


%Zoomed in plot:
figure(2)
plot(t,yFE)
hold on
plot(t,yBE)
plot(t,yRK4)
plot(t_ode45,y_ode45);
plot(t,ASval)
grid
legend('Forward Euler','Backward Euler','Runge Kutta 4th order', 'ode45','Analytical solution','Location','Northeast')
title('Simulation of world population using various methods')
xlabel('Year')
ylabel('population [million people]')
xlim([2002.68142,2002.68162])
ylim([6193.18,6193.205])
%xlim([2000,2010])

%  (f) error calculation
err_FE = yFE - ASval;
err_BE = yBE - ASval;
err_RK4 = yRK4 - ASval;
err_ode45 = y_ode45 - ASval';

fprintf('The maximum error of EF is %g\n',max(abs(err_FE)))
fprintf('The maximum error of BE is %g\n',max(abs(err_BE)))
fprintf('The maximum error of RK4 is %g\n',max(abs(err_RK4)))
fprintf('The maximum error of ode45 is %g\n',max(abs(err_ode45)))

figure(4)
hold on
% plot(t,err_FE)
% plot(t,err_BE)
%plot(t,err_RK4)
plot(t_ode45,err_ode45)
% legend('Forward Euler','Backward Euler','Runge Kutta 4th order', 'ode45')
% title('Error of the 4 methods')
% xlabel('Year')
% ylabel('Error [Million people]')


%Actual error:
% act_err_FE = yFE(1:10) - ASval;
% act_err_BE = yBE - ASval;
% act_err_RK4 = yRK4 - ASval;
% act_err_ode45 = y_ode45 - ASval_ode45;
% 
% fprintf('The actual maximum error of FE is %g\n',max(abs(act_err_FE)))
% fprintf('The actual maximum error of BE is %g\n',max(abs(act_err_BE)))
% fprintf('The actual maximum error of RK4 is %g\n',max(abs(act_err_RK4)))
% fprintf('The actual maximum error of ode45 is %g\n',max(abs(act_err_ode45)))



%% Problem 2
clc; clear variables; close all;




% (a)
% Explaining how to set up system and solving system using Runge-Kutta integrator


% (b)
% Actually solving for 0 <= t <= 25.

% IV's:
x0 = 1; y0 = 1; z0 = 1;
IC = [x0; y0; z0];

%time vector
h = 0.001; % time step
t = 0:h:25; % time vector

%Solving
Y = odesRK4(@ode_Lorenz, IC, t);
x_val = Y(:,1);
y_val = Y(:,2);
z_val = Y(:,3);

% PLOTTING
%xy, xz, and yz:
subplot(3, 2, 1); plot(x_val, y_val), xlabel('x'); ylabel('y'); title('xy-plane')
subplot(3, 2, 2); plot(x_val, z_val), xlabel('x'); ylabel('z'); title('xz-plane')
subplot(3, 2, 3); plot(y_val, z_val), xlabel('y'); ylabel('z'); title('yz-plane')


% x, y ,z versus t
subplot(3, 2, 4); plot(t, x_val); xlabel('t'); ylabel('x'); title('x vs. t')
subplot(3, 2, 5); plot(t, y_val); xlabel('t'); ylabel('y'); title('y vs. t')
subplot(3, 2, 6); plot(t, z_val); xlabel('t'); ylabel('z'); title('z vs. t')

%Finding endpoint:
h = 0.001; % time step
t2 = 0:h:100; % time vector

%Solving
Y = odesRK4(@ode_Lorenz, IC, t2);
x_end = Y(end,1);
y_end = Y(end,2);
z_end = Y(end,3);

fprintf('the endpoint which the solution is attracted to is (%g, %g, %g)',x_end,y_end,z_end)


% (c) 
% Changing r, and plotting again
%Solving
Y28 = odesRK4(@ode_Lorenz2, IC, t);
x_val28 = Y28(:,1);
y_val28 = Y28(:,2);
z_val28 = Y28(:,3);

figure(3)
% PLOTTING
%xy, xz, and yz:
subplot(3, 2, 1); plot(x_val28, y_val28); xlabel('x'); ylabel('y'); title('xy-plane')
subplot(3, 2, 2); plot(x_val28, z_val28); xlabel('x'); ylabel('z'); title('xz-plane')
subplot(3, 2, 3); plot(y_val28, z_val28); xlabel('y'); ylabel('z'); title('yz-plane')

% x, y ,z versus t
subplot(3, 2, 4); plot(t, x_val28); xlabel('t'); ylabel('x'); title('x vs. t')
subplot(3, 2, 5); plot(t, y_val28); xlabel('t'); ylabel('y'); title('y vs. t')
subplot(3, 2, 6); plot(t, z_val28); xlabel('t'); ylabel('z'); title('z vs. t')

% (d) 
% unpredictability for r = 28:

%To different initial points with IV's:
% IV's:
x02 = 6; y02 = 6; z02 = 6;
IC2 = [x02; y02; z02];

x03 = 6; y03 = 6.01; z03 = 6;
IC3 = [x03; y03; z03];

%Solving both:
Y_IC2 = odesRK4(@ode_Lorenz2, IC2, t);
x_IC2 = Y_IC2(:,1);
y_IC2 = Y_IC2(:,2);
z_IC2 = Y_IC2(:,3);

Y_IC3 = odesRK4(@ode_Lorenz2, IC3, t);
x_IC3 = Y_IC3(:,1);
y_IC3 = Y_IC3(:,2);
z_IC3 = Y_IC3(:,3);

figure (4)
 
plot3(x_IC2,y_IC2,z_IC2);
grid
hold on
plot3(x_IC3,y_IC3,z_IC3,'r');
xlabel('x')
ylabel('y')
zlabel('z')
legend('IV = (6,6,6)','IV=(6,6.01,6)')

figure(5)
subplot(1, 2, 1);
plot(-y_IC2,z_IC2);
legend('IV = (6,6,6)')
xlabel('y')
ylabel('z')
subplot(1, 2, 2);
plot(-y_IC3,z_IC3,'r');
xlabel('y')
ylabel('z')
legend('IV=(6,6.01,6)')
function dYdt = ode_Lorenz(Y,t)
% Return derivatives of Y at t
% Y = [x; y, z]
% dYdt = [dy/dt; dyD/dt,dyDD/dt]

% Given values:
sigma = 10;
b = 8/3;
r = 20;

x = Y(1,1);
y = Y(2,1);
z = Y(3,1);

xD = sigma*(y - x);
yD = r*x - y - x*z;
zD = x*y - b*z;

% build COLUMN vector of outputs
dYdt(1,1) = xD;
dYdt(2,1) = yD;
dYdt(3,1) = zD;

end
function dYdt = ode_Lorenz2(Y,t)
% Return derivatives of Y at t
% Y = [x; y, z]
% dYdt = [dy/dt; dyD/dt,dyDD/dt]

% Given values:
sigma = 10;
b = 8/3;
r = 28;

x = Y(1,1);
y = Y(2,1);
z = Y(3,1);

xD = sigma*(y - x);
yD = r*x - y - x*z;
zD = x*y - b*z;

% build COLUMN vector of outputs
dYdt(1,1) = xD;
dYdt(2,1) = yD;
dYdt(3,1) = zD;

end
