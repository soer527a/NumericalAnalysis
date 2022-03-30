function y = odeFE(dydt,y0,t)
% Solves an ordinary differential equation using the % Forward Euler method. %
% Inputs: dydt function handle for evaluating the ODE
%y0 initial condition/value
% t time vector [t_0 t_1 t_2 ... t_n]
% Output: y resulting values at each time step in t

%creating a y-vector
y = zeros(1,length(t));
y(1) = y0;

%calculating next values in y
for i = 1:(length(t)-1)
    %defining step length h (may vary)
    h = t(i+1)-t(i);

    y(1+i) = y(i) + h.*dydt(y(i),t(i));
end

end
