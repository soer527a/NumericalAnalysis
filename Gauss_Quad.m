function [int_Gauss,varargout] = Gauss_Quad(f,a,b, varargin)
% Gauss Quadrature integration calculator
% Computes the integral of f from a to b
% Uses the Gauss Quadrature rule to compute the integral
% f = function
% a = boundary 1
% b = boundary 2
% n = polynomial

%h is the width of each trapezoid.

if nargin == 4
    n = varargin{1};
    int_Gauss = Gauss_Quad_Calc(f,a,b,n);

elseif nargin == 5
    if strcmp(varargin{1},'Accuracy')
        A = varargin{2};

        Int_Exact = integral(f,a,b);
        for k_quad = 1:5
            Int_quad = Gauss_Quad_Calc(f,a,b,k_quad);
            E = Int_quad - Int_Exact;
            if  abs(E) < 10^(-A)
                break
            end

        end
        int_Gauss = Int_quad;
        varargout{1} = k_quad;

    end
end

    function int_Gauss_Calc = Gauss_Quad_Calc(f,a,b,n)
        J = (b-a)/2;
        [xi,w] = GaussTable(n);

        x = (b+a)/2 + ((b-a)/2) * xi;

        result = 0;
        for i = 1:n
            result = result + w(i)*f(x(i));
        end
        int_Gauss_Calc = result * J;
    end

end