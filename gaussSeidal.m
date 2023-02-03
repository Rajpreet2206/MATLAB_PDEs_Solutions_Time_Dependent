%% Gauss-Seidel Solver
% The discretized equation for the Heat Equation is given by:
% [T(i+1,j) -2*T(i,j) + T(i-1,j)]/hx^2 + [T(i,j+1) -2*T(i,j) + T(i,j-1)]/hy^2 = f(xi,yi)
% where i,j denote the x and y nodes respectively
% 
% This can be simplified to:
% T(i,j) = [ hx^2*hy^2 / 2*(hx^2 + hy^2) ] * { [T(i+1,j) + T(i-1,j)]/hx^2 + [T(i,j+1) + T(i,j-1)]/hy^2 - f(xi,yi) }

%% INPUTS: 
% b: Vector of KNOWN variables
% Nx: Number of x nodes with unknown temperatures
% Ny: Number of y nodes with unknown temperatures
% tol: Tolerance for Gauss-Seidal Solver

%% OUTPUTS:
% T_noB: Matrix with temperatures at nodes excluding boundaries
% storage: size rquirement of this method
%%
function [T_noB, storage] = gaussSeidal(T_curr,Nx,Ny,tol,P)
    % Constants
    
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    k = (hx^2 * hy^2) / (2 * (hx^2 + hy^2));
    k = (1/(1 + 2*dt/hx^2 + 2*dt/hy^2))
    % Initialize Temperature Matrix && Residual
    
    T_next = zeros(Ny+2,Nx+2);
    % To assign non-zero boundaries
    %  T(first row) = y_max boundary
    %  T(last row) = y_min boundary values
    %  T(first_column) = x_min boundary vlues
    %  T(last_column) = x_max boundary values
    
    A_x = zeros(size(T));
    b_k = zeros(size(T));
    res_norm = 1;
    
    % Convergence Loop
    while res_norm > tol
        
        % Nodal Loop
        for i=2:(Nx+1)
            for j=2:(Ny+1)
                x = (i-1)*hx;
                y = (j-1)*hy;
                % Gauss-Seidal Method
                T_next(j,i) = k * ( dt*( T_next(i-1,j) + T_next(i+1,j) )/hy^2 + dt*( T_next(i,j-1) + T_next(i,j+1) )/hx^2 + T_curr(j,i));
                
                % For Residual Calculations
                b_k(j,i) = b(x,y);
            end
        end
        
        % Check for convergence
        for i=2:(Nx+1)
            for j=2:(Ny+1)                
                A_x(j,i) = T_next(i,j( T_next(i-1,j) + T_next(i+1,j) -2*T(i,j) )/hy^2 + ( T(i,j-1) + T(i,j+1) -2*T(i,j) )/hx^2;
            end
        end
        
        res_norm = sqrt(sum(sum((b_k - A_x).^2))/(Nx*Ny));
  
        % For next iteration
        Told = T;
        % For Output with No Boundary Nodes
        T_noB = T(2:Ny+1,2:Nx+1);
    end
    
%% This section is only implemented once, and not implemented for the Timeit function.
    
    if(P==1)
        %Calculating Storage Requirement
        
        storage = numel(b)+numel(Nx)+numel(Ny)+numel(tol)+numel(1)+numel(hx)+numel(hy)+...
                  +numel(k)+numel(T)+numel(Told)+numel(res_norm)+numel(i)+numel(j)+numel(x)+numel(y)+...
                  +numel(A_x)+numel(b_k)+numel(T_noB);
              
        % Creating plots
        title = strcat('Gauss Seidel Solver for Nx = Ny = ', num2str(Nx));

        plotter(T,Nx,Ny,title);
    end 

end
