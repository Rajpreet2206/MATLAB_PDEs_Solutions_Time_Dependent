function [T_next] = explicit_euler(Nx,Ny,dt,T_curr)

    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    N = Nx*Ny;
    A = zeros(N:N); % system matrix
    T_next = zeros(Ny,Nx); % Matrix of node temperatures

    %[A,b_known_terms] = matrixA(Nx,Ny);
    A = matrixA_exp(Nx,Ny);

    T_next = T_curr + dt*(A*T_curr);
end