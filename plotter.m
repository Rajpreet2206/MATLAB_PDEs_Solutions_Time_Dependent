% INPUTS
% T: Temperature values, including boundaries

function [f] = plotter(T,Nx,Ny,what)
    x = linspace(0,1,Nx+2);
    y = linspace(0,1,Ny+2);


    % Create surface plot

    f = figure('Visible','off');
    surf(x,y,T);
    xlabel("X axis"); ylabel("Y axis");
    filename = sprintf("%s.png",what);
    %fullname = fullfile('.\plots', filename);
    saveas(f,filename);
end