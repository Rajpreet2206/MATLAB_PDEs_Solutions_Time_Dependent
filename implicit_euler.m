function [T_next] = implicit_euler(Nx,Ny,dt,T_curr)

    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    N = Nx*Ny;
    T_next = zeros(Ny,Nx); % Matrix of node temperatures
    T_next = GaussSeidel(Nx,Ny,dt,T_curr);
    end