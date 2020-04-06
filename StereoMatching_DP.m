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

dsi = zeros(height, height);
dsiSize = size(dsi,2);
disparity = zeros(height, width);

for x=1:width
    tic;
    im_near_subs = zeros(window, window, height);
    im_far_subs = zeros(window, window, height);
    
    parfor idx=1:height
        im_near_subs(:, :, idx) = imgLookAt(im2double(im_near_360), LON(idx, x), LAT(idx, x), window, fov);
        im_far_subs(:, :, idx) = imgLookAt(im2double(im_far_360), LON(idx, x), LAT(idx, x), window, fov);
    end

    for nearY=1:height
        im_near_sub = im_near_subs(:, :, nearY);
        for farY=1:height
            im_far_sub = im_far_subs(:, :, farY);
            diff = double(im_near_sub) - double(im_far_sub);
            diffSq = diff .^2;
            SSE = sum(sum(diffSq));
            dsi(farY, nearY) = SSE;
        end
    end
    dsi = dsi./max(max(dsi));
    
    C = zeros(dsiSize, dsiSize);
    Pointers = zeros(dsiSize, dsiSize);
    occlusionConstant = 0.1;
    
    for i = 2 : dsiSize
        C(i, 1) = C(i-1,1) + occlusionConstant;
        C(1, i) = C(1,i-1) + occlusionConstant;
        Pointers(i,1) = 2;
        Pointers(1,i) = 3;
    end
    
    for i = 2 : dsiSize
        for j = 2 : dsiSize
            values = [C(i-1,j-1) + dsi(i,j), C(i-1,j) + occlusionConstant, C(i,j-1) + occlusionConstant];
            [val, idx] = min(values);
            C(i, j) = val;
            Pointers(i,j) = idx;
        end
    end
    i = dsiSize;
    j = dsiSize;
    path = [];
    
    while(i ~= 1 || j ~= 1)
        switch(Pointers(i,j))
            case 1
                i = i -1;
                j = j - 1;
            case 2
                i = i - 1;
            case 3
                j = j - 1;
        end
        path = [path;i,j];
    end
    
    newScanline = zeros(1, dsiSize);
    
    for zX = 1 : size(path,1)
        xLoc = path(zX,1);
        yLoc = path(zX,2);
        if(zX < size(path,1) && xLoc == path(zX+1,1) + 1 && yLoc == path(zX+1,2) + 1)
        newScanline(xLoc) = sqrt((xLoc - yLoc) ^ 2);      
        else
        newScanline(xLoc) = 0; 
        end
    end
    
    disparity(:, x) = newScanline;
    toc;
    disp(['width : ' num2str(x)]);
end
disparity = disparity./max(max(disparity));



