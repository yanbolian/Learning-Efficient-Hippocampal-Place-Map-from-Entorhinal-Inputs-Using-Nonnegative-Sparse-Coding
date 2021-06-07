function array = display_matrix(A, resize_factor, num_column, figure_on)
%  This function display the columns of A in an image where each rectangle
%  of the image represents the column of A
% A: may represent the connections between two populations;
% resizeFactor: the factor for a higher resolution
% numColumns: # of rectangles in each row of the image
% onFlag: whether to display the image

[L0, M] = size(A); % L: length of the column; M: number of columns of A
sz0 = ceil( sqrt(L0) ); % The smallest size of a square large enough to display each column

% Zero-padded version of A
A_Padded = zeros(sz0^2, M);
A_Padded(1:L0, 1:M) = A;
% min(A_Padded(:))

% Set defaul value for 'resizeFactor'
if ~exist('resize_factor','var')
    resize_factor=1;
end

sz = sz0 * resize_factor; % size of the resized square
L = sz^2; % 
A_Padded = reshape(imresize(reshape(A_Padded,sz0,sz0,M),resize_factor),L,M);
A_Padded(A_Padded<0) = 0;
% min(A_Padded(:))

% By default, the displayed image is square unless defined by 'numColumns'
if ~exist('num_column', 'var')
    nCol = ceil(sqrt(M));
else
    nCol = num_column;
end
nRow = ceil(M/nCol); % number of rows of the displayed image

buff = 2 * resize_factor; % thickness of the boundary between different blocks in the image
black = 1.02; % Color 

% The array that saves the pixel values of the displayed image
array = black * ones( buff + nRow * (sz+buff), buff + nCol * (sz+buff) ); 

% fill 'array ' column by column
m=1; % index of the column of A
for i = 1 : nRow
    for j = 1 : nCol
        if m < M + 1
            % clim: larget absolute value in m-th column of A
            clim = max( abs( A_Padded(:,m) ) );
            if clim == 0
                clim = 1;
            end
            % fill m-th column of A into 'array'
            array(  buff + (i-1) * (sz+buff) + (1:sz), ...
                    buff + (j-1) * (sz+buff) + (1:sz) ) = ...
                    reshape(A_Padded(:,m), sz, sz) / clim;
            
            % used when plotting samples
%             array(  buf + (i-1) * (sz+buf) + (1:sz), ...
%                     buf + (j-1) * (sz+buf) + (1:sz) ) = ...
%                     reshape(A_Padded(:,m), sz, sz);
        end
        m = m + 1; % update m with the index of next column
    end
end

% The image will be displayed by default unless stated otherwise
if ~exist('figure_on', 'var')
    figure_on = 1;
end

% Display the image
if figure_on == 1
%     colormap gray
    imagesc(array,[min(min(array(:)),0), black]);
    axis image off
    
%     surf(array);
%     shading interp
%     view(0,180)
%     colormap gray
% %     axis off
%     axis image
%     drawnow
end

end
