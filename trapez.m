function [Int_T,varargout] = trapez(f,a,b,varargin)
% Trapezoidal integration calculator
% Computes the integral of f from a to b using n+1 points
% Uses the composite trapezoidal rule to compute the integral
% Added possibility of various types of plots with varargin: 'Trap_Plot' or
% 'Anim Plot'. With only 4 inputs, computes only the integral

% f = function
% a = boundary 1
% b = boundary 2
% n = number of panels

%h is the width of each trapezoid.

%---------------- NORMAL CALCULATION WITH SPECIFIED n ---------------
if nargin == 4
    n = varargin{1};
    [h,Int_T] = trapez_calc(f,a,b,n);
end

% ------------- DIFFERENT POSSIBILITIES WITH 5 INPUTS -------------
if nargin == 5
    if strcmp(varargin{1},'Accuracy')
        A = varargin{2} + 1;
    else
        n = varargin{1};
         [h,Int_T] = trapez_calc(f,a,b,n);
    end

    % -- INTEGRAL CALCULATION WITH SPECIFIED NUMBER OF ACCURATE DECIMALS --
    if strcmp(varargin{1},'Accuracy')

        Int_Exact = integral(f,a,b);
        for k_trap = 1:10^A
             [h, Int_Trap] = trapez_calc(f,a,b,k_trap);
            E = Int_Trap - Int_Exact;
            if  abs(E) < 10^(-A)
                break
            end

        end
        Int_T = Int_Trap;
        varargout{1} = k_trap;





        % ----- PLOTTING INTEGRAL CALCULATION WITH DIFFERENT OPTIONS ------
        %CHOOSING A NORMAL PLOT WITH TRAPEZOIDS
    elseif strcmp(varargin{2},'Trap_Plot')

        %Plotting the Patches
        for j = 0:(n-1)
            x_val = [a + j*h, a + j*h, a + h + j*h, a + h + j*h];
            y_val = [0, f(a + j*h), f(a + h + j*h), 0];

            patch(x_val,y_val,'green','HandleVisibility','off')
            hold on
        end

        %Plotting the Function
        fplot(f,'k','LineWidth',1.5)
        xlim([a,b])
        title('Trapezoidal Integral')
        legend('function f')

        %CHOOSING ANIMATED PLOT WITH AREA OF EACH PANEL (NOT RECOMMENDED FOR HIGH
        %n-VALUES)
    elseif strcmp(varargin{2},'Anim_Plot')

        %Plotting the Function
        fplot(f,'k','LineWidth',1.5)
        xlim([a,b])
        title('Trapezoidal Integral')
        legend('function f')
        hold on

        %Plotting the Patches and integral values
        for j = 0:(n-1)
            pause(0.5)
            x_val = [a + j*h, a + j*h, a + h + j*h, a + h + j*h];
            y_val = [0, f(a + j*h), f(a + h + j*h), 0];

            patch(x_val,y_val,'green','HandleVisibility','off')

            panel_area = (h/2) * (f(a + j*h) + f(a + h + j*h));
            area_text = strcat('Area =',string(panel_area));
            text(a+(j+0.5)*h,f(a+(j+0.5)*h),area_text)




        end

        Total_Area_txt = strcat('Integral = ',string(Int_T));
        text(a + (1/2) * (b-a) , f(a + (b-a)/10),Total_Area_txt,'FontSize',15,'FontName','Times New Roman')


    end
end


% ------------- GENERAL TRAPEZOIDAL CALCULATION FUNCTION ------------
    function [h,Integral] = trapez_calc(f,a,b,n)
        h = (b-a)/n;
        result2 = 0.5*f(a) + 0.5*f(b);
        for L = 1:(n-1)
            result2 = result2 + f(a + L*h);
        end
        Integral = h*result2;
    end


end