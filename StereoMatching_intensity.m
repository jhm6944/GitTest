addpath('PanoBasic-master');
add_path;

im_near = rgb2gray(imread('IMG_A0.jpg'));
im_far   = rgb2gray(imread('IMG_A1.jpg'));

im_near = imresize(im_near, [size(im_near, 1)/4, size(im_near, 2)/4]);
im_far = imresize(im_far, [size(im_far, 1)/4, size(im_far, 2)/4]);

im_near_360 = zeros(640, 1280, 'uint8');
im_far_360 = zeros(640, 1280, 'uint8');

im_near_360(:, 481:800) = im_near;
im_far_360(:, 481:800) = im_far;

[TX, TY] = meshgrid(1:size(im_near_360, 2), 1:size(im_near_360, 1));
LON_TEMP = (TX/(size(im_near_360, 2)/2) - 1) .* pi();
LAT_TEMP = -1*(TY/(size(im_near_360, 1)/2) - 1) .* pi()/2;

LON = LON_TEMP(:, 481:800);
LAT = LAT_TEMP(:, 481:800);

width = size(im_near, 2);
height = size(im_near, 1);
window = 16;
fov = pi()/12;

featureSize = (window ^ 2);
disparity = zeros(height, width);
cost = zeros(height, width);
for x=1:width
    tic;
    nearDescriptors = zeros(featureSize, height);
    farDescriptors = zeros(featureSize, height);
    disparityLine = zeros(height, 1);
    costLine = zeros(height, 1);
    
    for y=1:height
        im_near_sub = imgLookAt(im2double(im_near_360), LON(y, x), LAT(y, x), window, fov);
        nearDescriptors(:, y) = im_near_sub(:);
        im_far_sub = imgLookAt(im2double(im_far_360), LON(y, x), LAT(y, x), window, fov);
        farDescriptors(:, y) = im_far_sub(:);
        
%         subplot(2, 1, 1);
%         imshow(im_near_sub);
%         subplot(2, 1, 2);
%         imshow(im_far_sub);
%         pause;
    end
    
    bound = 40;
    for farY=1:height
        SSDs = zeros(bound, 1);
        count = 0;
        if(farY < height/2)
            for nearY = farY:-1:max(farY-bound, 1)
                count = count + 1;
                diff = nearDescriptors(:, nearY) - farDescriptors(:, farY);
                diffSq = diff .^2;
                SSDs(count, 1) = sum(sum(diffSq));
            end 
            SSDs_temp = SSDs(1:count, 1);
            [val, loc] = min(SSDs_temp);
            disparityLine(farY, 1) = loc;
            costLine(farY, 1) = val;
        elseif(farY > height/2)
            for nearY = farY:min(farY+bound, height)
                count = count + 1;
                diff = nearDescriptors(:, nearY) - farDescriptors(:, farY);
                diffSq = diff .^2;
                SSDs(count, 1) = sum(sum(diffSq));
            end
            SSDs_temp = SSDs(1:count, 1);
            [val, loc] = min(SSDs_temp);
            disparityLine(farY, 1) = loc;
            costLine(farY, 1) = val;
        else
            disparityLine(farY, 1) = 0;
            costLine(farY, 1) = 0;
        end
    end
    
    disparity(:, x) = disparityLine;
    cost(:, x) = costLine;
    toc;
    disp(['width : ' num2str(x)]);
end
disparity = disparity./max(max(disparity));


