function color_map = jet_modified(num)

if ~ exist('num', 'var')
    num = 256;
end

k = 0; % background level

c = jet(num);

color_map = c;
background_range = floor(k*num);

color_map(1:num, :) = c;
color_map(1:background_range, :) = repmat(c(1,:),background_range,1);
color_map(num+1, :) = [0 0 0]; % black for the boundaries