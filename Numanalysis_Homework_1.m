%% Numerical Analysis Homework 1

%% Problem 1 - Newton Raphson
clc, clear variables, close all;
format long
fprintf('<strong>Results of Problem 1 </strong> \n\n')

%using polar coordinates to plot the circle. we know r = sqrt(5):
r = sqrt(5);
theta = linspace(0,2*pi,100);

x_circle = cos(theta)*r;

y_circle = sin(theta)*r;

figure(1)
plot(x_circle,y_circle)
hold on

x_parabola = linspace(-3,3,100);
y_parabola = x_parabola.^2 - 1; 


plot(x_parabola,y_parabola)
ylim([-3,3])
axis equal

%Graphically estimated intersections are x  = -1.60606 and x = 1.60606

%To use Newton-Raphson i rewrite the two equations. we get:

% y_parabola - y_circle = 0.

f = @(x) -sqrt(5-x.^2)+(x.^2-1);
%fp = @(x) 2*x+(x./(sqrt(5-x.^2)));
p = @(x) x.^2 - 1;
n=0;
for x = [-1,1]

while abs(f(x))>eps  &&  n<20
%x = x - f(x)/fp(x);
h = 0.001;
    fpn = (f(x+h)-f(x))/h;
    x = x - f(x)/fpn;

n = n + 1;

plot(x,p(x),'r*')
hold on
end
fprintf('an intersection was found to be %g after %g iterations\n',x,n)
n = 0;
end


%% Problem 2 - data fitting
clc;
clear variables; close all

fprintf('<strong>Results of Problem 2 </strong> \n\n')


%Defining the full matrices of years and toxin levels.
Year = [1993, 1995, 1997, 1999, 2001, 2003, 2005, 2007];
Tox = [12.0, 12.7, 13.0, 15.2,  18.2, 19.8, 24.1, 28.1];
xx = linspace(1993 ,2009,1001);

%creating lagrange interpolation data for on these data values
yint = Lagrange(Year,Tox,xx);

%Plotting the interpolation
figure(2)
plot(Year,Tox,'ko',"MarkerFaceColor","r")
xlabel('Year')
ylabel('Toxin conc.')
xlim([1993, 2009])
ylim([0,max(yint)*1.1])
hold on
plot(xx,yint,'b-','LineWidth',1.5)

title('Lagrange Interpolation')

%Printing resulting prediction for 2009
fprintf("the predicted toxin concentration in 2009 is %g \n",Lagrange(Year,Tox,2009))

%Defining new matrix with data, leaving out 1997 and 1999
Year2 = [1993, 1995, 2001, 2003, 2005, 2007];
Tox2 = [12.0, 12.7,  18.2, 19.8, 24.1, 28.1];

xx2 = linspace(1993 ,2007,1001);
yint2 = Lagrange(Year2,Tox2,xx2);

figure(3)
plot(Year2,Tox2,'ko',"MarkerFaceColor","r")
xlabel('x')
ylabel('y')
xlim([1993, 2007])
ylim([min(yint2)*1.1,max(yint2)*1.1])
hold on
plot(xx2,yint2,'b-','LineWidth',1.5)
title('Lagrange Interpolation')

%Printing predictions from 1997 and 1999
fprintf("\nthe predicted toxin concentration in 1997  is %g \n",Lagrange(Year2,Tox2,1997))
fprintf("the predicted toxin concentration in 1999  is %g \n",Lagrange(Year2,Tox2,1999))


%Repeating b) with spline instead
%creating lagrange interpolation data for these data values
yspline = spline(Year,Tox,xx);

%Plotting the interpolation
figure(4)
plot(Year,Tox,'ko',"MarkerFaceColor","r")
xlabel('Year')
ylabel('Toxin conc.')
xlim([1993, 2009])
ylim([0,max(yspline)*1.1])
hold on
plot(xx,yspline,'b-','LineWidth',1.5)

title('Cubic spline Interpolation')

%Printing resulting prediction for 2009
fprintf("the predicted toxin concentration in 2009 is %g \n",spline(Year,Tox,2009))


%% Problem 3 - Integral calculations.

clc;
clear variables; %close all;
fprintf('<strong>Results of Problem 3 </strong> \n\n')

% Task:
% Evaluate the given integral with 4 different methods:
% - Analytically
% - Trapezoidal rule
% - Simpson's rule
% - Gauss Quadrature

% The function to be integrated:
f = @(x) -0.0547.*x.^4 + 0.8646.*x.^3 - 4.1562.*x.^2 + 6.2917.*x + 2;

%Boundaries a and b of integration:
a = 0;
b = 8;

%Exact integral
Exact_int = integral(f,a,b);
fprintf('The exact integral is                  %g\n',Exact_int) 

%Integration with trapezoidal rule½
Int_T = trapez(f,a,b,8);
fprintf('The integral using trapezoidal rule is %g\n',Int_T)

%Integration with simpson
Int_Simpson = simpson(f,a,b,8);
fprintf('The integral using simpson rule is     %g\n',Int_Simpson)

%Integration with Gauss Quadrature
Int_Gauss = Gauss_Quad(f,a,b,3);
fprintf('The integral using gauss quadrature is %g\n',Int_Gauss)

%Accurate trapezoidal
[Int_Trap, n_Trap] = trapez(f,a,b,'Accuracy',3);
fprintf('\nThe accurate Integral using trapezoidal is %g with n = %g\n', Int_Trap, n_Trap)

%Accurate simpson
[Int_SimpA, n_Simp] = simpson(f,a,b,'Accuracy',3);
fprintf('The accurate Integral using simpson is %g with n = %g\n', Int_SimpA, n_Simp)

%Accurate Gauss Quadrature
[Int_GaussA, n_Gauss] = Gauss_Quad(f,a,b,'Accuracy',3);
fprintf('The accurate Integral using Gauss Quadrature is %g with n = %g\n', Int_GaussA, n_Gauss)