function y = odeBE(dydt, y0, t)
% solves an ordinary differential equation using the Backward Euler method.
%
% Inputs: dydt(y,t)  Function handle for evaluating the ODE
%         y0         Initial condition
%         t          Time vector [t_0 t_1 t_2 ... t_n]
% Output: y          Resulting function values at each time step in t

y = zeros(1,length(t));
y(1)=y0;
dy = 10^-6;

for i = 1:(length(t)-1)

    % current time step length (may vary)
    h = t(i+1)-t(i);

    y_next = y(i); % initial guess for next point
    g = 10; %high function value to get loop started
   


    while abs(g)>0.00001  &&  i<20
         g = y_next - (y(i) + h * dydt(y_next,t(i+1)));
         g2 = y_next + dy - (y(i) + h * dydt(y_next + dy, t(i+1)));
         gp = (g2 - g)/dy;

         %newton-raphson:
         y_next = y_next  - g/gp;

    end
y(i+1) = y_next;
end

end