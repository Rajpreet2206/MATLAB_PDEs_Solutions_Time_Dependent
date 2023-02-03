function [Tn]=GaussSeidel(Nx,Ny,dt,T_curr)

    tol = 1.0e-06;
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    k = (1 + 2*dt/hx^2 + 2*dt/hy^2);
    N = Nx*Ny;
    Tp = T_curr;
    Tn = zeros(Ny+2,Nx+2);
    A_x = zeros(size(Tn));
    res_norm = 1.0;
   
    while res_norm > tol
     
     % Nodal Loop
        for i=2:(Nx+1)
            for j=2:(Ny+1)
                % Gauss-Seidal Method
                Tn(j,i) = 1/k * ( dt*( Tp(i-1,j) + Tn(i+1,j) )/hy^2 + dt*( Tp(i,j-1) + Tn(i,j+1) )/hx^2 + T_curr(j,i));
            end
        end
       
         % Check for convergence
        for i=2:(Nx+1)
            for j=2:(Ny+1)                           
                A_x(j,i) =  Tn(i,j)*k - dt*( Tn(i-1,j) + Tn(i+1,j) )/hy^2 - dt*( Tn(i,j-1) + Tn(i,j+1) )/hx^2;            
            end
        end
        
        % Residual Norm Computation
        res_norm = sqrt(sum(sum((T_curr - A_x).^2))/(Nx*Ny));
  
        % For next iteration
        Tp = Tn;
    end
end