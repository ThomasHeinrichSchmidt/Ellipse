% Try to use VLFeat


close all;
clear all;
debug_on_warning(1);
debug_on_error(1);


% img = im2double ( imread ('http://t-h-schmidt.de/Filipp/stimulus2.jpg'));
img = im2double(imread('stimulus2.jpg'));
img(img>0)=1;

% stop for debugging
% keyboard()

% Define meshgrid
[destX,destY] = meshgrid(1:size(img,2),1:size(img,1));

(x1, x2) := grid(I);
phi := vl_tps(I);
(xp1, xp2) := vl_wtps(phi, annRef);
for x := 1 to length(size(I, 1)) do
for y := 1 to length(size(I, 2)) do
    I1(x, y) = I(xp1(x, y), xp2(x, y));
end
end


% Plot ellipse, bend ellipse, and sine wave ellipse
close all;

figure
subplot(2,2,1)
show(img)
axis xy
colormap gray
title 'Original'

subplot(2,2,2)
show(I1)
axis xy
colormap gray
title 'transformed'

