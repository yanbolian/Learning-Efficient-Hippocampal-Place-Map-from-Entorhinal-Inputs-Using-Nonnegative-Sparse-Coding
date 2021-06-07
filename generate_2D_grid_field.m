function G = generate_2D_grid_field(x_size, y_size, Nx, Ny, lambda, orien, phase, figure_on)
% Generate a 2D grid field using the sum of three sine gratings (Solstad et al. 2006)
% 
% G: returned 2D grid field normalised to [0, 1] with size Ny * Nx (Ny rows and Nx columns)
% x_size, y_size: size of the environment (unit: m)
% Nx, Ny: number of discrete points in the x and y axis
% lambda: grid spacing (unit: m)
% orien: orientation of the grid pattern (unit: radians); corresponding to the vertical line in a clockwise manner
% phase: grid phase - (x0,y0) (unit: m)

freq = 2/(sqrt(3)*lambda); % corresponding frequency; derived from lambda
x0 = phase(1);
y0 = phase(2);

% Generate XY coordinates in the environment
[X, Y] = meshgrid(linspace(0,x_size,Nx), linspace(0,y_size, Ny)); 

% Grid fields with values normalized to [0, 1]
G = 1/3 *   (...
        2/3 *   (...
            cos(2*pi*freq * ((X-x0)*cos(2/3*pi*1+orien) + (Y-y0)*sin(2/3*pi*1+orien)))...
            + cos(2*pi*freq * ((X-x0)*cos(2/3*pi*2+orien) + (Y-y0)*sin(2/3*pi*2+orien)))...
            + cos(2*pi*freq * ((X-x0)*cos(2/3*pi*3+orien) + (Y-y0)*sin(2/3*pi*3+orien)))...
                )...
        + 1 ...
            );

if exist('figure_on','var')
    if figure_on==1
    figure;
    surf(X, Y, G);
    shading interp
    view(0,90)
    colormap jet
%     axis off
    axis image
    end
end