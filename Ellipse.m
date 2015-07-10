% As of Filipp's email 27.6.2014 
% Die % bezeichnen Kommentare.
% meshgrid und interp2 brauche ich, um Transformationen anwenden zu können.
% bwboundaries, cell2mat und poly2mask brauche ich, um Bilder in Vektoren umzuwandeln und umgekehrt.

% Uses geom2D toolbox
% Uses contourinfoplot function (http://ruccs.rutgers.edu/~manish/demos/curveinfo.html)
% Uses sinefunction function

% uses package image for: im2double()
pkg load image
% may also load all packages
%pkg load all

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

sourceX = destX;
sourceY = destY;

img_trans = interp2(destX,destY,img,sourceX,sourceY,'cubic');
img_trans(isnan(img_trans)==1)= 0;

% Convert image to boundary vectors
img_trans_sine = bwboundaries(img_trans);
img_trans_sine = cell2mat(img_trans_sine);

dummy(:,1) = img_trans_sine(:,2); dummy(:,2) = img_trans_sine(:,1);
img_trans_sine(:,1)=dummy(:,1);img_trans_sine(:,2)=dummy(:,2);

% Find maximas of curvature by using function contourinfoplot
[surprisalmap] = contourinfoplot(img_trans_sine(1:end,1)', img_trans_sine(1:end,2)', 0, 5);
% Find points of maximal curvature by searching within MANUALLY defined areas of the shape !!!!!!!!!!!!!!!!!
[C1,I1]=max(surprisalmap(1:100,1));    % search in maxima map between position 1:100
max1 = img_trans_sine(I1,:);
[C2,I2]=max(surprisalmap(400:450,1));  % search in maxima map between position 300:600
max2 = img_trans_sine(I2+300,:);

% Define part of img_trans_sine between maxima max1 and max2 as "roi"
roi=img_trans_sine(I1:I2+300,:);

% Define sine wave
x1 = [0:1:max(roi(:,1))-min(roi(:,1))];               % adjust sine wave to lenght of roi
oscill_s = 5 + 0.5;                                     % number of oscillations + 0.5 oscillations
freq_s = oscill_s * 2*pi / (max(roi(:,1))-min(roi(:,1))); % frequency
amp_s = 20;                                             % amplitude
phase_s = 0;                                            % phase
offset_s = 0;                                           % offset
[ z ] = sinefunction( x1, amp_s, freq_s, phase_s, offset_s ); % define sine wave function
x2 = [min(roi(:,1)):1:max(roi(:,1))];       % translocate sine wave to x values of roi
plot(x1,z)
sinevector(:,1)=x2;
sinevector(:,2)=z';

% Define triangular wave with the same characteristics as the sine wave

st = (max(roi(:,1))-min(roi(:,1)))/(oscill_s*4);

steps_dummy = [+amp_s -amp_s -amp_s  +amp_s;   % steps on x axis
               +st    +st    +st     +st  ];   % steps on y axis

steps = zeros(2,oscill_s*4-0.5*4);
for j=1:(oscill_s-0.5)+1
    steps(:,(j*4-3)+1:(j*4)+1)=steps_dummy;
end
steps_dummy = [-2*amp_s ; 0];
steps = [steps steps(:,end)+steps_dummy(:,1)]; % steps(:,end)+steps_dummy(:,2)];

triavector_dummy = zeros(2,oscill_s*4);
for i=1:(oscill_s*4)+1
    triavector_dummy(2,i+1)=triavector_dummy(2,i)+steps(1,i);
    triavector_dummy(1,i+1)=triavector_dummy(1,i)+steps(2,i);
end

triavector_dummy=triavector_dummy';
triavector_dummy(1,:)=[];
plot(triavector_dummy(:,1),triavector_dummy(:,2))

y_output=[];
for xxx = 1:oscill_s*4
    for xx = round(triavector_dummy(xxx,1)):round(triavector_dummy(xxx+1,1))
y_dummy(xx+1)=triavector_dummy(xxx,2)+((triavector_dummy(xxx+1,2)-triavector_dummy(xxx,2))/(triavector_dummy(xxx+1,1)-triavector_dummy(xxx,1))) *(xx-triavector_dummy(xxx,1));
    end
end

triavector(:,1)=x2';
triavector(:,2)=y_dummy';

% Calculate new y values of a morphvector between sinevector and triavector

weight_sine = 0;   % set relative weights of sine function 0...1 so that weight_sine + weight_tria = 1
weight_tria = 1;   % set relative weights of triangular function 0...1 so that  weight_sine + weight_tria = 1

morphvector(:,1) = x2';
morphvector(:,2) = [( weight_sine*(sinevector(:,2)+amp_s) + weight_tria*(triavector(:,2)+amp_s) )/(weight_sine+weight_tria)];

% Calculate the normal vector of neighbours for every point of roi

% Example for roi(2,:)
% dx = [roi(3,1)-roi(1,1)];
% dy = [roi(3,2)-roi(1,2)];
% normal = [roi(2,1)-dy roi(2,2)+dx];

n = 20; % define the distance to the neighbours ("smooth")

roi_new=zeros(length(roi),2);

for k=1:(length(roi)-n*2)  % no calculation of normal vectors for first and last point
    k=k+n;

    % linear function defined by the two neighbours
    m_neigh = (roi(k+n,2)-roi(k-n,2))/(roi(k+n,1)-roi(k-n,1));
    if m_neigh ==Inf || m_neigh ==-Inf
        m_neigh=0
    end
    b_neigh = -1*m_neigh*roi(k-n,1)+roi(k-n,2);
    % linear function of normal vector passing through roi(k,:)
    m = -1/m_neigh;
    if m ==Inf || m ==-Inf
        m=0
    end
    b = -1*m*roi(k,1)+roi(k,2);

    % Calculate the intersections of the normal vector function and a circle
    % with center [v w] = [ roi(k,1) roi(k,2) ] and the radius
    % r = abs(morphvector(find(morphvector==(roi(k,1))),2))
    v = roi(k,1);
    w = roi(k,2);
    r = abs(morphvector(find(morphvector==(roi(k,1))),2));

    vprim = (1 + m^2);
    wprim = 2 * m * (b - w) - 2 * v;
    uprim = v^2 + (b - w)^2 - r^2;

    delta = wprim^2 - 4 * vprim * uprim

    x1_intersect = (-wprim + sqrt(delta)) / (2 * vprim);
    y1_intersect = m * x1_intersect + b;
    intersect1 = [x1_intersect y1_intersect];

    x2_intersect = (-wprim - sqrt(delta)) / (2 * vprim);
    y2_intersect = m * x2_intersect + b;
    intersect2 = [x2_intersect y2_intersect];

    % Test whether the resulting intersections x1,y1 and x2,y2 are on the
    % left (or on the right) of the line between k-1 and k+1 by calculating
    % the cross product A x B = Ax * By - Ay * Bx
    cross_product = @(u,o) ((u(1,1)*o(1,2))-(u(1,2)*o(1,1)));

    point_left_of_line = cross_product(roi(k+n,:)-roi(k-n,:)+roi(k,:), intersect1(1,:)-roi(k,:) )

    % Assign intersect1 or intersect2 depending on result of cross product
    if morphvector(find(morphvector==(roi(k,1))),2) > 0 && point_left_of_line > 0
        roi_new(k,:) = intersect1(1,:);
    elseif morphvector(find(morphvector==(roi(k,1))),2) < 0 && point_left_of_line > 0
        roi_new(k,:) = intersect2(1,:);
    elseif morphvector(find(morphvector==(roi(k,1))),2) > 0 && point_left_of_line < 0
        roi_new(k,:) = intersect2(1,:);
    elseif morphvector(find(morphvector==(roi(k,1))),2) < 0 && point_left_of_line < 0
        roi_new(k,:) = intersect1(1,:);
    end

    % Assign y values of sinevector in straight parts of the contour
    if m == 0
        roi_new(k,1) = roi(k,1);
        roi_new(k,2) = roi_new(k,2)+morphvector(find(morphvector==(roi(k,1))),2);
    end

end

% Define start and end points that do not have any normal vectors
for kk=1:n
    roi_new(kk,1) = roi(kk,1);
    roi_new(kk,2) = roi(kk,2);
    roi_new(end-kk+1,1) = roi(end-kk+1,1)
    roi_new(end-kk+1,2) = roi(end-kk+1,2)
end

% Insert transformed roi into img_trans_sine
img_trans_sine(I1:I2+300,:)=roi_new(:,:);
img_trans_sine = poly2mask(img_trans_sine(:,1)',img_trans_sine(:,2)',600,600);

% Plot ellipse, bend ellipse, and sine wave ellipse
close all;

figure
subplot(2,2,1)
show(img)
axis xy
colormap gray
title 'original'

subplot(2,2,2)
show(img_trans)
axis xy
colormap gray
title 'transformed'

subplot(2,2,3)
show(img_trans_sine)
axis xy
colormap gray
title 'transformed' 

