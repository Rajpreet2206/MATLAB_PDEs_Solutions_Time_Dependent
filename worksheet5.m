%% Worksheet 4: Solving Stationary Partial Differential Equations

%% General Instructions  
% All the plots will be displayed after all the grid sizes are completed
%% Defining Constants and Functions
clearvars;
close all;
clear all;
clc;

rmdir plots s;
mkdir plots;

%% Solving Full Matrix, Sparse Matrix and Gauss Seidel for various Grid Sizes

iter = 0;
dt_array = [1/64,1/128,1/256,1/512,1/1024,1/2048,1/4096];
N_array = [3,7,15,31];
ts_array = [1/8,2/8,3/8,4/8];

ii = 1;

for Grid_size = N_array
   for dt = dt_array
     t = 0.0;
     disp('Calculating results by Explicit Method for grid size below: Please Wait !');
     disp(Grid_size);
     
     iter = iter+1;
     Nx = Grid_size; Ny = Grid_size;
     N = Nx*Ny;
     x = linspace(0,1,Nx+2);
     y = linspace(0,1,Ny+2);
     T0 = ones(N,1);
     T_inner = zeros(Ny,Nx);
     T_curr = T0;
     kk = 1; 
      while (t < 4/8)
          T_exp = explicit_euler(Nx,Ny,dt,T_curr);
          T_curr = T_exp;
          t = t + dt;
          for tsnap = ts_array
              if (t == tsnap)
                  for i = 1 : Nx
                    for j = 1 : Ny
                        iter = (i-1)*Ny + j;
                        T_inner(j,i) = T_exp(iter,1);
                    end
                  end
                  T_bound = [zeros(1,Nx); T_inner; zeros(1,Nx)];
                  T_bound = [zeros(Ny+2,1), T_bound, zeros(Ny+2,1)];
           
                  %format rat;
                  %plottitle = sprintf('Surface Plot: Nx=Ny=%d dt=%s at t=%s',Nx,strtrim(rats(dt)),strtrim(rats(tsnap)));
                  plottitle = sprintf('Surface Plot: Nx=Ny=%d dt=%s at t=%s',Nx,dt,tsnap);
                  plotter(T_bound,Nx,Ny,plottitle);
                  figure(kk);
                  subplot(4,7,ii);
                  surf(x,y,T_bound);
                  kk = kk + 1;
                  break;
              end

          end

      end

      ii = ii+1;

   end
  
end

    for i = 1:length(ts_array)
                    figure(i)
                    annotation('textbox',[0.45,0.005,0.1,0.1],'String',"Surface Plots at t="+rats(ts_array(i)))
                    annotation('textbox',[0.05,0.75,0.1,0.1],'String',"Nx=Ny=" + N_array(1))
                    annotation('textbox',[0.05,0.55,0.1,0.1],'String',"Nx=Ny=" + N_array(2))
                    annotation('textbox',[0.05,0.35,0.1,0.1],'String',"Nx=Ny=" + N_array(3))
                    annotation('textbox',[0.05,0.12,0.1,0.1],'String',"Nx=Ny=" + N_array(4))
                    annotation('textbox',[0.15,0.88,0.1,0.1],'String',"dt="+rats(dt_array(1)))
                    annotation('textbox',[0.26,0.88,0.1,0.1],'String',"dt="+rats(dt_array(2)))
                    annotation('textbox',[0.37,0.88,0.1,0.1],'String',"dt="+rats(dt_array(3)))
                    annotation('textbox',[0.48,0.88,0.1,0.1],'String',"dt="+rats(dt_array(4)))
                    annotation('textbox',[0.59,0.88,0.1,0.1],'String',"dt="+rats(dt_array(5)))
                    annotation('textbox',[0.70,0.88,0.1,0.1],'String',"dt="+rats(dt_array(6)))
                    annotation('textbox',[0.81,0.88,0.1,0.1],'String',"dt="+rats(dt_array(7)))
    end

%Implicit Method
kimp_ini = kk;
ii = 1;
for Grid_size = N_array
    t = 0.0;
    dt = 1/64;
    Nx = Grid_size; Ny = Grid_size;
    N = Nx*Ny;
    T0 = ones(Ny,Nx);
    TB0 = [zeros(1,Nx); T0; zeros(1,Nx)];
    TB0 = [zeros(Ny+2,1), TB0, zeros(Ny+2,1)];
    T_curr = TB0;
    x = linspace(0,1,Nx+2);
    y = linspace(0,1,Ny+2);
    kk = kimp_ini;


    while (t < ts_array(end) )
          T_imp = implicit_euler(Nx,Ny,dt,T_curr);
          T_curr = T_imp;
          t = t + dt;
      for tsnap = ts_array
        if (t == tsnap)
                 
                  figure(kk);
                  %format rat;
                  %figtitle = sprintf('Implicit Method (Gauss Seidel):Surface Plots for dt %f at t=%f',rats(dt),rats(tsnap));
                  figtitle = sprintf('Implicit Method (Gauss Seidel):Surface Plots for dt %f at t=%f',dt,tsnap);

                  sgtitle(figtitle);

                  subplot(2,2,ii);
                  surf(x,y,T_imp);
                  xlabel("X axis"); ylabel("Y axis");
                  plottitle = sprintf('Nx = Ny = %d',Nx);
                  title(plottitle);
                  kk = kk + 1;
                  break;
        end
      end
    end
    ii = ii + 1;
end
    
 






