function [Int_Simpson, varargout] = simpson(f,a,b, varargin)
% Simpson rule integration calculator
% Computes the integral of f from a to b using n+1 points
% NUMBER OF PANELS n MUST BE EVEN
% Uses the composite trapezoidal rule to compute the integral
% f = function
% a = boundary 1
% b = boundary 2
% n = number of panels


%---------------- NORMAL CALCULATION WITH SPECIFIED n ---------------
if nargin == 4
    n = varargin{1};
    if n/2 == floor(n/2)
        % if n is even the integral is calculated
        
        [h,Int_Simpson] = Simpson_Calc(f,a,b,n);

    else
        %if n is odd the integral is not calculated
        fprintf('error, n is not an even number')
        Int_Simpson = 'Not Calculated';
    end
end

if nargin == 5
    if strcmp(varargin{1},'Accuracy')
        A = varargin{2} + 1;

        % --- calculating the integral for different n's -----
        Int_Exact = integral(f,a,b);
        for k_simp = 1:10^A
             [h, Integral] = Simpson_Calc(f,a,b,k_simp);
            E = Integral- Int_Exact;
            if  abs(E) < 10^(-A)
                break
            end

        end
        Int_Simpson = Integral;
        varargout{1} = k_simp;


    end
end

%---------- GENERAL SIMPSON RULE CALCULATION FUNCTION ------------
    function [h,integral] = Simpson_Calc(f,a,b,n)
        %
        %h is the width of each interval
        h = (b-a)/n;
        result = (f(a)+f(b));
        for i = 1:2:(n-1)
            result = result + 4*f(a + i*h);
        end
        for j = 2:2:(n-2)
            result = result + 2 * f(a + j*h);
        end
        integral = (h/3)*result;
    end

end

