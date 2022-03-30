function y = odeRK4(dydt, y0, t)
% Solves an ordinary differential equation using the
% 4th order Runge-Kutta method.
%
% Inputs: dydt function handle for evaluating the ODE
% y0 initial condition/value
% t time vector [t_0 t_1 t_2 ... t_n]
% Output: y resulting values at each time step in t

%creating a y-vector
y = zeros(1,length(t));
y(1)=y0;


%calculating next values in y
for i = 1:(length(t)-1)
 
    %defining step length h (may vary)
    h = t(i+1)-t(i);

k1 = dydt(y(i), t(i));
k2 = dydt(y(i)+(h/2)*k1,t(i) + h/2 );
k3 = dydt(y(i)+(h/2)*k2,t(i) + h/2 );
k4 = dydt(y(i) + h*k3, t(i) + h);

y(1+i) = y(i) + (h/6)*k1 + (h/3)*(k2 + k3) + (h/6) * k4;

end

end