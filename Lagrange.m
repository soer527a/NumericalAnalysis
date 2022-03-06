function yint = Lagrange(x,y,xx)
% Lagrange interpolating polynomial
% yint= Lagrange(x,y,xx): Uses an (n -1)-order
% Lagrange interpolating polynomial based on n data points
% to determine a value of the dependent variable (yint) at
% a given value of the independent variable, xx.
% inputs:
% x = independent variable at known points
% y = dependent variable at known points
% xx = value(s) of independent variable at which the interpolation is done
% output:
% yint= interpolated value(s) of dependent variable


n = length(x);
n_int = length(xx);
yint=zeros(n_int,1);

for k = 1:n_int
    Y = 0;
    for j = 1:n
        L = 1;
        for i = 1:n
            if i~=j
                L = L * ((xx(k) - x(i))/(x(j)-x(i)));
            end
        end
        Y = Y + y(j) * L;
    end
    yint(k)= Y;
end