function [xi,w] = GaussTable(nGP)
% Returns sampling points and weights for for Gauss quadrature of order n
%
% INPUTS:
%   nGP:      order of quadrature
%
% OUTPUTS:
%    xi:      sampling points
%     W:      weight functions


    switch nGP 
        case 1
            xi = 0;
            w  = 2;
        case 2
            xi = [-sqrt(1/3) sqrt(1/3)];
            w  = [1 1];
        case 3
            xi = [-sqrt(3/5) 0 sqrt(3/5)];
            w  = [5/9 8/9 5/9];
        case 4
            xi(1) = -sqrt(3/7 +2/7*sqrt(6/5));
            xi(2) = -sqrt(3/7 -2/7*sqrt(6/5));
            xi(3) =  sqrt(3/7 -2/7*sqrt(6/5));
            xi(4) =  sqrt(3/7 +2/7*sqrt(6/5));
            w(1) = (18-sqrt(30))/36;
            w(2) = (18+sqrt(30))/36;
            w(3) = (18+sqrt(30))/36;
            w(4) = (18-sqrt(30))/36;
        case 5
            xi(1) = -1/3*sqrt(5+2*sqrt(10/7));
            xi(2) = -1/3*sqrt(5-2*sqrt(10/7));
            xi(3) = 0;
            xi(4) =  1/3*sqrt(5-2*sqrt(10/7));
            xi(5) =  1/3*sqrt(5+2*sqrt(10/7));
            w(1) = (322-13*sqrt(70))/900;
            w(2) = (322+13*sqrt(70))/900;
            w(3) = 128/225;
            w(4) = (322+13*sqrt(70))/900;
            w(5) = (322-13*sqrt(70))/900;
    end

end