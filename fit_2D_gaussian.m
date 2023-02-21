%% Fit to Gaussian field

function [fields_fit, fit_parameters] = fit_2D_gaussian(A, x_size, y_size, Nx, Ny, num_peak)


opts.resize_factor = 3; % the upsampling factor for the bases prior to fitting
opts.Nx = Nx * opts.resize_factor; % image patches size in pixels
opts.Ny = Ny * opts.resize_factor; % image patches size in pixels
opts.num_cells = size(A, 2); % Number of cells to be fitted
opts.plot_figure = 1; % if set to 1, then for each filter the original and fitted data will be plotted
opts.num_trials = 500; % number of trials in optimization loop
opts.fit_error_threshold = 5; % normalized fitting error threshold in percentage
opts.datetime = datetime; % Date and time of the fitting process

x = linspace(0,x_size,opts.Nx);
y = linspace(0,y_size,opts.Ny);
[X, Y] = meshgrid(x, y);
XY(:, :, 1) = X;
XY(:, :, 2) = Y;

if ~exist('num_peak', 'var')
    num_peak = 1;
end

% p = a * exp{ -ln5 [(x-x0)^2+(y-y0)^2)] / r^2 }; 1/5 response at radius
% params(1) = amplitude; % Amplitude of the 2D Gaussian
% params(2,3) = x0, y0; % Center of the 2D Gaussian
% params(4) = radius; % Radius of 2D Gaussian
fun_2D_gaussian = @(params, XY) ...
    params(1) * exp( -log(5) * ((XY(:, :, 1)-params(2)).^2+(XY(:, :, 2)-params(3)).^2) / params(4)^2 );

for i_cell = 1 : opts.num_cells
    
    fprintf('Cell #%d. Starting 2D Gaussian fitting\n', i_cell);
    
    field = reshape(A(:, i_cell), Nx, Ny); % convert each basis from vector into matrix
	field_resized = imresize(field, opts.resize_factor, 'lanczos3'); % upsampling
    params = [];
    fit_error = 100;
    
    for i_trials = 1 : opts.num_trials
        fprintf('(%d)', i_trials);
        params_initial = rand(4*num_peak, 1);
        options = optimset('Display', 'off', 'MaxIter', 1000);
        
        if max(field_resized(:)) == 0
            params = [0 0 0 1e-15];
            break
        end
        
        [params_current, resnorm] = lsqcurvefit(fun_2D_gaussian, ...
            params_initial, XY, field_resized, [], [], options);
        fit_error_current = resnorm / sum(field_resized(:).^2) * 100; % compute the fit error percentage
        
        % update fitError and determine whether to exit fitting process
        if fit_error_current <= fit_error
            fit_error = fit_error_current;
            params = params_current;
        end
        if fit_error < opts.fit_error_threshold
            break;
        end
        
        % print
        if ~mod(i_trials, 10)
            fprintf(' min fit error reached = %.2f percent.\n', fit_error);
        end
    end
    fprintf('Done!\n');
%     params
    params(4) = abs(params(4));
    field_fit = fun_2D_gaussian(params, XY);
    
    if opts.plot_figure
        figure(1);
        clf
        subplot(121);
        maxi = max(abs(field_resized(:))) + 1e-15;
        imagesc(x, y, field_resized, [-maxi maxi]); title('Place field')
        colormap(gray)
        colorbar
        xlabel('m'); ylabel('m');
        axis image
        
        drawnow
        subplot(122);
        maxi = max(abs(field_fit(:))) + + 1e-15;
        imagesc(x, y, field_fit, [-maxi maxi]); title('Fitted Filter')
        colormap(gray)
        colorbar
        xlabel('m'); ylabel('m');
        axis image
        
        drawnow
    end

    % save parameters
    fields_fit(i_cell).params = params;
    fields_fit(i_cell).gaussian_fitted = field_fit;
    fields_fit(i_cell).fit_error = fit_error;
    
    params'
    fprintf('The normalized fitting error is %.2f percent.\n', fit_error);
    
%     pause
end
fit_parameters = opts;