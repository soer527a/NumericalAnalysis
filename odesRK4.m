function Y = odesRK4(dYdt,Y0,t)
% Runge-Kutta 4th order method for solving a system of ODE's.
%
% INPUTS: dYdt(Y,t) = function returning the derivatives (column vector)
% Y0 = Initial conditions (column vector)
% t = Time [t_0; t_1; t_2; ...; t_n]
%
% OUTPUTS: Y = Matrix with resulting function values over time
% Y(n,i) = i'th function value at n'th timestep, i.e. each
% row contains all function values for that time

% get size of problem, no. of equations, timesteps. etc
ne = length(Y0); % no. equations
nt = length(t);  % timesteps
h  = t(2)-t(1);  % step size


% initialize output
Y = zeros(nt,ne);
Y(1,:) = Y0(:)'; % store results row-wise

%defining step length h
h = max(t)/(length(t)-1);

%calculating next values in Y, for each timestep t
for n = 1:(length(t)-1)
    Yn = Y(n,:)';
    tn = t(n);

    k1 = dYdt(Yn, tn);
    k2 = dYdt(Yn+(h/2)*k1,tn + h/2 );
    k3 = dYdt(Yn+(h/2)*k2,tn + h/2 );
    k4 = dYdt(Yn+h*k3, tn + h);

    Ynext = Yn + (h/6)*k1 + (h/3)*(k2 + k3) + (h/6) * k4;
    Y(n+1, :) = Ynext';

end
end


