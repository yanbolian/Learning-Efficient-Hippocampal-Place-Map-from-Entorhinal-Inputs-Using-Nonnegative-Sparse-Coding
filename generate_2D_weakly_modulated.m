function G = generate_2D_weakly_modulated(x_size, y_size, Nx, Ny, figure_on)
% Generate a 2D grid field using the sum of three sine gratings (Solstad et al. 2006)
% 
% G: returned 2D grid field normalised to [0, 1] with size Ny * Nx (Ny rows and Nx columns)
% x_size, y_size: size of the environment (unit: m)
% Nx, Ny: number of discrete points in the x and y axis

G0 = rand(Nx, Ny);
sigma_analog = 0.06 + 0.01*randn; % unit: m; the smoothing Gaussian kernel in Neher et al. 2017
sigma_discrete = sigma_analog / x_size * Nx; % Assume x_size:y_size = Nx:Ny

G = imgaussfilt(G0, sigma_discrete);

G = (G - min(G(:))) / (max(G(:)) - min(G(:))) ;

if exist('figure_on','var')
    if figure_on==1
    figure;
    imagesc(G);
    colormap(jet)
    colorbar
    axis image
    end
end