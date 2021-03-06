%% Set Params
addpath('image-morphing-master');
addpath('PanoBasic-master');
add_path;

global PARAMS;

PARAMS.LFU_W = 50;
PARAMS.HEIGHT = 2560;
PARAMS.WIDTH  = 5120;
PARAMS.TIMES = 500;

PARAMS.RESIZE_HEIGHT = PARAMS.HEIGHT/4;
PARAMS.RESIZE_WIDTH = PARAMS.WIDTH/4;

[TX, TY] = meshgrid(1:PARAMS.RESIZE_WIDTH, 1:PARAMS.RESIZE_HEIGHT);
TX = TX(:);
TY = TY(:);

LF = cell(2, 1);
LF_FILE = cell(2, 1);
% LF_FILE{1, 1} = ['Row7'; 'Row8'; 'Row9'];
% LF_FILE{2, 1} = ['Row13'; 'Row14'; 'Row15'];

LF_FILE{1, 1} = ['Row15'; 'Row16'; 'Row17'];
LF_FILE{2, 1} = ['Row21'; 'Row22'; 'Row23'];

for i=1:length(LF)
    LF{i, 1} = zeros(3*PARAMS.LFU_W * PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3, 1, 'uint8');
    for n=1:PARAMS.LFU_W*3
        if(n <= PARAMS.LFU_W)
            file = sprintf('../Matlab_FeatureExtraction/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(1, :), n);
        elseif(n <= PARAMS.LFU_W*2)
            file = sprintf('../Matlab_FeatureExtraction/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(2, :), n - PARAMS.LFU_W);
        else
            file = sprintf('../Matlab_FeatureExtraction/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(3, :), n - PARAMS.LFU_W*2);
        end
        img = imread(file);
        img = imresize(img, [PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH]);
        img_C1 = img(:, :, 1);
        img_C2 = img(:, :, 2);
        img_C3 = img(:, :, 3);
        
%         if(i == 2)
%             nn = (PARAMS.LFU_W * 3 - n);
%         else
%             nn = n - 1;
%         end
        nn = n - 1;
        
        LF{i, 1}((nn) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF{i, 1}((nn) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF{i, 1}((nn) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 1) = img_C3(:);
        
        disp(['LF ' num2str(i) ' - ' num2str(n) ', ' num2str(nn)]);
    end
end
%% DEBUG FEATURES

height = PARAMS.HEIGHT/4;
width = PARAMS.WIDTH/4;
windowsize = 50;
fov = pi()/12;

LFU_W = PARAMS.LFU_W * 3;
ANGLE_S = -1 * deg2rad(45);
ANGLE_E = deg2rad(45);

h_t = ((height*2)-1) / floor(height*2/10);
hh = 1:h_t:(height*2);
 
OUT_W = floor((ANGLE_E - ANGLE_S) / deg2rad(360/width));
w_t = (OUT_W-1) / floor(OUT_W/10);
ww = 1:w_t:OUT_W;

for idx_w=1:size(ww, 2)
    w = ww(1, idx_w);
    angle = ANGLE_S + deg2rad(360/width) .* (w-1);
    P_near = ((PARAMS.LFU_W / 2) - 0) .* tan(angle) + 0 + (LFU_W / 2);
    P_far = -1 * ((PARAMS.LFU_W / 2) + 0) .* tan(angle) + 0 + (LFU_W / 2);
    P_near = round(P_near);
    P_far = round(P_far);
    
    U = angle * (180.0 / pi)*(1.0 / 180.0)*width / 2 + width / 2;
    U = round(U);
    
    im_near = rgb2gray(LoadFrame(LF{1, 1}, P_near));
    im_far  = rgb2gray(LoadFrame(LF{2, 1}, P_far));

    for idx_h=1:size(hh, 2)
        if(hh(1, idx_h) < height)
            idx_h_near = idx_h - (disparity(idx_h, idx_w) - 1);
        elseif(hh(1, idx_h) > height)
            idx_h_near = idx_h + (disparity(idx_h, idx_w) - 1);
        end
        h_near = hh(1, idx_h_near) - height/2;
        h_far = hh(1, idx_h) - height/2;
        
        w_near = U;
        r_near = 0;
        if(h_near < 1)
            h_near = -1 * h_near;
            w_near = U + width/2;
            r_near = 1;
        elseif(h_near > height)
            h_near = 2 * height - h_near;
            w_near = U + width/2;
            r_near = 1;
        end
        if(w_near > width)
            w_near = w_near - width;
        end
        
        w_far = U;
        r_far = 0;
        if(h_far < 1)
            h_far = -1 * h_far;
            w_far = U + width/2;
            r_far = 1;
        elseif(h_far > height)
            h_far = 2*height - h_far;
            w_far = U + width/2;
            r_far = 1;
        end
        if(w_far > width)
            w_far = w_far - width;
        end
        
        lon_near = (w_near/(width/2) - 1) * pi();
        lat_near = -1*(h_near/(height/2) - 1) * pi()/2;
        lon_far  = (w_far/(width/2) - 1) * pi();
        lat_far  = -1*(h_far/(height/2) - 1) * pi()/2;
        
        im_near_sub = imgLookAt(im2double(im_near), lon_near, lat_near, windowsize, fov);
        im_near_sub(isnan(im_near_sub)) = 0;
        if(r_near == 1)
            im_near_sub = rot90(im_near_sub, 2);
        end
        im_far_sub = imgLookAt(im2double(im_far), lon_far, lat_far, windowsize, fov);
        im_far_sub(isnan(im_far_sub)) = 0;
        if(r_far == 1)
            im_far_sub = rot90(im_far_sub, 2);
        end
        
        disp([num2str(idx_h_near) ', ' num2str(idx_h) ', ' ...
            num2str(w_near) ', ' num2str(h_near) ', ' num2str(w_far) ', ' ...
            num2str(h_far)]);
        
        subplot(1, 2, 1);
        imshow(im_near_sub);
        subplot(1, 2, 2);
        imshow(im_far_sub);
        pause;
    end
end

%% Run (w blending)
[posx, posy] = meshgrid(-24:12:24, -24:12:24);
posx = posx(:);
posy = posy(:);

% for idx = 1:size(posx, 1)
%     curr_posx = posx(idx, 1); 
%     curr_posy = posy(idx, 1);
    curr_posx = 0;
    curr_posy = 0;
    posx_front = curr_posx;            posy_front = curr_posy;

    morph_image = renderingLF_morph(LF{1, 1}, LF{2, 1}, posx_front, posy_front);
    imshow(morph_image);
% end

%%
function MORPHED_IMAGE = renderingLF_morph(LF_A0, LF_A1, POS_X, POS_Y)%, MatchedFeatures)
    global PARAMS;
    LFU_W = PARAMS.LFU_W * 3;
    ANGLE_S = -1 * deg2rad(45);
    ANGLE_E = deg2rad(45);
    OUT_W = floor((ANGLE_E - ANGLE_S) / deg2rad(360/PARAMS.RESIZE_WIDTH));
    VIEW_W = (1:OUT_W)';
    ANGLE = ANGLE_S + deg2rad(360/PARAMS.RESIZE_WIDTH) .* (VIEW_W-1);
        
    P_near = ((PARAMS.LFU_W / 2) - POS_Y) .* tan(ANGLE) + POS_X + (LFU_W / 2);
    P_far = -1 * ((PARAMS.LFU_W / 2) + POS_Y) .* tan(ANGLE) + POS_X + (LFU_W / 2);

    P_near = round(P_near);
    P_far = round(P_far);
    P_near(P_near < 1) = 1;    P_near(P_near > LFU_W) = LFU_W;
    P_far(P_near < 1) = 1;    P_far(P_near > LFU_W) = LFU_W;
    
    U_FRONT = ANGLE.*(180.0 / pi)*(1.0 / 180.0)*PARAMS.RESIZE_WIDTH / 2 + PARAMS.RESIZE_WIDTH / 2;
    U_BACK = U_FRONT + PARAMS.RESIZE_WIDTH/2;
    
    U_FRONT = round(U_FRONT);
    U_BACK = round(U_BACK);
    U_FRONT(U_FRONT > PARAMS.RESIZE_WIDTH) = U_FRONT(U_FRONT > PARAMS.RESIZE_WIDTH) - PARAMS.RESIZE_WIDTH;
    U_BACK(U_BACK > PARAMS.RESIZE_WIDTH) = U_BACK(U_BACK > PARAMS.RESIZE_WIDTH) - PARAMS.RESIZE_WIDTH;
    
    H = (1:PARAMS.RESIZE_HEIGHT)';
    
    P_near = repmat(P_near', size(H, 1), 1);
    P_far = repmat(P_far', size(H, 1), 1);
    U_FRONT = repmat(U_FRONT', size(H, 1), 1);
    U_BACK = repmat(U_BACK', size(H, 1), 1);
    
    VIEW_FRONT_NEAR_3 = im2uint8(inter8_mat(LF_A0, 0, P_near-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 1));
    VIEW_FRONT_NEAR_2 = im2uint8(inter8_mat(LF_A0, 0, P_near-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 2));
    VIEW_FRONT_NEAR_1 = im2uint8(inter8_mat(LF_A0, 0, P_near-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 3));
    VIEW_FRONT_NEAR = cat(3, VIEW_FRONT_NEAR_1, VIEW_FRONT_NEAR_2, VIEW_FRONT_NEAR_3);
    
    VIEW_BACK_NEAR_3 = im2uint8(inter8_mat(LF_A0, 0, P_near-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 1));
    VIEW_BACK_NEAR_2 = im2uint8(inter8_mat(LF_A0, 0, P_near-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 2));
    VIEW_BACK_NEAR_1 = im2uint8(inter8_mat(LF_A0, 0, P_near-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 3));
    VIEW_BACK_NEAR = cat(3, VIEW_BACK_NEAR_1, VIEW_BACK_NEAR_2, VIEW_BACK_NEAR_3);
    
%     view_near = cat(1, rot90(VIEW_BACK_NEAR(1:size(VIEW_BACK_NEAR, 1)/2, :, :), 2), ...
%         VIEW_FRONT_NEAR, rot90(VIEW_BACK_NEAR((size(VIEW_BACK_NEAR, 1)/2 + 1):size(VIEW_BACK_NEAR, 1), :, :), 2));

    view_near = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
    view_near(:, (PARAMS.RESIZE_WIDTH/2 - PARAMS.RESIZE_WIDTH/8 + 1):(PARAMS.RESIZE_WIDTH/2 + PARAMS.RESIZE_WIDTH/8), :) = VIEW_FRONT_NEAR;
    view_near(:, 1:PARAMS.RESIZE_WIDTH/8, :) = VIEW_BACK_NEAR(:, (size(VIEW_BACK_NEAR, 2)/2 + 1):size(VIEW_BACK_NEAR, 2), :);
    view_near(:, (PARAMS.RESIZE_WIDTH - PARAMS.RESIZE_WIDTH/8 + 1):PARAMS.RESIZE_WIDTH, :) = VIEW_BACK_NEAR(:, 1:size(VIEW_BACK_NEAR, 2)/2, :);
    
    VIEW_FRONT_FAR_3 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 1));
    VIEW_FRONT_FAR_2 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 2));
    VIEW_FRONT_FAR_1 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 3));
    VIEW_FRONT_FAR = cat(3, VIEW_FRONT_FAR_1, VIEW_FRONT_FAR_2, VIEW_FRONT_FAR_3);
    
    VIEW_BACK_FAR_3 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 1));
    VIEW_BACK_FAR_2 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 2));
    VIEW_BACK_FAR_1 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 3));
    VIEW_BACK_FAR = cat(3, VIEW_BACK_FAR_1, VIEW_BACK_FAR_2, VIEW_BACK_FAR_3);
    
%     view_far = cat(1, rot90(VIEW_BACK_FAR(1:size(VIEW_BACK_FAR, 1)/2, :, :), 2), ...
%         VIEW_FRONT_FAR, rot90(VIEW_BACK_FAR((size(VIEW_BACK_FAR, 1)/2 + 1):size(VIEW_BACK_FAR, 1), :, :), 2));

    view_far = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
    view_far(:, (PARAMS.RESIZE_WIDTH/2 - PARAMS.RESIZE_WIDTH/8 + 1):(PARAMS.RESIZE_WIDTH/2 + PARAMS.RESIZE_WIDTH/8), :) = VIEW_FRONT_FAR;
    view_far(:, 1:PARAMS.RESIZE_WIDTH/8, :) = VIEW_BACK_FAR(:, (size(VIEW_BACK_FAR, 2)/2 + 1):size(VIEW_BACK_FAR, 2), :);
    view_far(:, (PARAMS.RESIZE_WIDTH - PARAMS.RESIZE_WIDTH/8 + 1):PARAMS.RESIZE_WIDTH, :) = VIEW_BACK_FAR(:, 1:size(VIEW_BACK_FAR, 2)/2, :);

%     showMatchedFeatures(view_near, view_far, features(:, 1:2), features(:, 3:4));

    height = PARAMS.HEIGHT/4;
    width = PARAMS.WIDTH/4;
    
    windowsize = 50;
    fov_1 = pi()/12;
    fov_2 = pi()/12;
    
    OUT_W = floor((ANGLE_E - ANGLE_S) / deg2rad(360/width));
    h_t = ((height*2)-1) / floor(height*2/10);
    hh = 1:h_t:(height*2);
    w_t = (OUT_W-1) / floor(OUT_W/10);
    ww = 1:w_t:OUT_W;
    
    disparity = zeros(size(hh, 2), size(ww, 2));
    cost = zeros(size(hh, 2), size(ww, 2));
    
    for idx_w=1:size(ww, 2)
        disp(['idx_w : ' num2str(idx_w) '/' num2str(size(ww, 2))]);
        tic; 
        w = ww(1, idx_w);
        angle = ANGLE_S + deg2rad(360/width) .* (w-1);
        U = angle * (180.0 / pi)*(1.0 / 180.0)*width / 2 + width / 2;
        U = round(U);
        
        im_near = rgb2gray(view_near);
        im_far  = rgb2gray(view_far);
        
        im_near_subs_1 = zeros(windowsize, windowsize, size(hh, 2));
        im_far_subs_1  = zeros(windowsize, windowsize, size(hh, 2));
        im_near_subs_2 = zeros(windowsize, windowsize, size(hh, 2));
        im_far_subs_2  = zeros(windowsize, windowsize, size(hh, 2));
        
        disparityLine = zeros(size(hh, 2), 1);
        costLine = zeros(size(hh, 2), 1);
        for idx_h=1:size(hh, 2)
            h = hh(1, idx_h) - height/2;
            w_n = U;
            rot = 0;
            if(h < 1)
                h = -1 * h;
                w_n = U + width/2;
                rot = 1;
            elseif(h > height)
                h = 2 * height - h;
                w_n = U + width/2;
                rot = 1;
            end
            
            if(w_n > width)
                w_n = w_n - width;
            end
            
            lon = (w_n/(width/2) - 1) * pi();
            lat = -1*(h/(height/2) - 1) * pi()/2;

%             disp([num2str(lon) ', ' num2str(lat)]);
            im_near_sub = imgLookAt(im2double(im_near), lon, lat, windowsize, fov_1);
            im_near_sub(isnan(im_near_sub)) = 0;
            im_near_sub_2 = imgLookAt(im2double(im_near), lon, lat, windowsize, fov_2);
            im_near_sub_2(isnan(im_near_sub_2)) = 0;
            if(rot == 1)
                im_near_sub = rot90(im_near_sub, 2);
                im_near_sub_2 = rot90(im_near_sub_2, 2);
            end
            im_near_subs_1(:, :, idx_h) = im_near_sub;
            im_near_subs_2(:, :, idx_h) = im_near_sub_2;
            
            im_far_sub = imgLookAt(im2double(im_far), lon, lat, windowsize, fov_1);
            im_far_sub(isnan(im_far_sub)) = 0;
            im_far_sub_2 = imgLookAt(im2double(im_far), lon, lat, windowsize, fov_2);
            im_far_sub_2(isnan(im_far_sub_2)) = 0;
            if(rot == 1)
                im_far_sub = rot90(im_far_sub, 2);
                im_far_sub_2 = rot90(im_far_sub_2, 2);
            end
            im_far_subs_1(:, :, idx_h) = im_far_sub;
            im_far_subs_2(:, :, idx_h) = im_far_sub_2;
%             im_near_sub = im_near((h-50:h+50), (w_n-50):(w_n+50));
%             im_far_sub = im_far((h-50:h+50), (w_n-50):(w_n+50));
            
%             disp(num2str(idx_h));
%             subplot(1, 2, 1);
%             imshow(im_near_sub);
%             subplot(1, 2, 2);
%             imshow(im_far_sub);
%             pause;
        end

        bound = 10;
%         bound_sub = 2;
        for idx_h_far=1:size(hh, 2)
            farY = hh(1, idx_h_far);
            
            if(idx_h_far >= 90 && idx_h_far <= 104)
                disparityLine(idx_h_far, 1) = 1;
                costLine(idx_h_far, 1) = -1;
            elseif(farY < height)
                cost_buff = zeros(bound, 3);
                count = 0;
                for idx_h_near=idx_h_far:-1:max(idx_h_far - bound, 1)
                    count = count + 1;
                    [cost_i, cost_g, cost_f] = CalculateCost(idx_h_near, idx_h_far, im_near_subs_1, im_far_subs_1, windowsize);
                    cost_buff(count, :) = [cost_i cost_g cost_f];
                end
                cost_t = cost_buff(1:count, :);
                cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 250 * cost_t(:, 3);
                [val, loc] = min(cost_t2);
%                 cost_buff = zeros(bound, 3);
%                 count = 0;
%                 for idx_h_near = max(loc - bound_sub, 1):min(loc + bound_sub, bound)
%                     count = count + 1;
%                     [cost_i, cost_g, cost_f] = CalculateCost(idx_h_near, idx_h_far, im_near_subs_2, im_far_subs_2, windowsize);
%                     cost_buff(count, :) = [cost_i cost_g cost_f];
%                 end
%                 cost_t = cost_buff(1:count, :);
%                 cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 250 * cost_t(:, 3);
%                 [val, loc] = min(cost_t2);
                disparityLine(idx_h_far, 1) = loc;
                costLine(idx_h_far, 1) = val;
            elseif(farY > height)
                cost_buff = zeros(bound, 3);
                count = 0;
                for idx_h_near = idx_h_far:min(idx_h_far + bound, size(hh, 2))
                    count = count + 1;
                    [cost_i, cost_g, cost_f] = CalculateCost(idx_h_near, idx_h_far, im_near_subs_1, im_far_subs_1, windowsize);
                    cost_buff(count, :) = [cost_i cost_g cost_f];
                end
                cost_t = cost_buff(1:count, :);
                cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 250 * cost_t(:, 3);
                [val, loc] = min(cost_t2);
%                 cost_buff = zeros(bound, 3);
%                 count = 0;
%                 for idx_h_near = max(loc - bound_sub, 1):min(loc + bound_sub, bound)
%                     count = count + 1;
%                     [cost_i, cost_g, cost_f] = CalculateCost(idx_h_near, idx_h_far, im_near_subs_2, im_far_subs_2, windowsize);
%                     cost_buff(count, :) = [cost_i cost_g cost_f];   
%                 end
%                 cost_t = cost_buff(1:count, :);
%                 cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 250 * cost_t(:, 3);
%                 [val, loc] = min(cost_t2);
                disparityLine(idx_h_far, 1) = loc;
                costLine(idx_h_far, 1) = val;
            else
                disparityLine(idx_h_far, 1) = 1;
                costLine(idx_h_far, 1) = -1;
            end
        end
        toc;
        disparity(:, idx_w) = disparityLine;
        cost(:, idx_w) = costLine;
    end
    
%     load('disparity');
%     load('cost');
    features = zeros(0, 5);
    for x=1:size(disparity, 2)    
        features_temp = zeros(0, 5);
        for y=1:size(disparity, 1)
            w = ww(1, x);
            h = hh(1, y);
            
            if(h < height)
                h_n = hh(1, y - (disparity(y, x) - 1));
            elseif(h > height)
                h_n = hh(1, y + (disparity(y, x) - 1));
            else
                h_n = h;
            end

            features_temp = cat(1, features_temp, [w, h_n, w, h, cost(y, x)]);
        end
        [C,ia,ic] = unique(features_temp(:, 1:2), 'rows');
        features = cat(1, features, features_temp(ia, :));
    end
    






    ratio = (POS_Y + 50)/100;
    
    pts_mean = (features(:, 1:2) + features(:, 3:4)) / 2;
    [C,ia,ic] = unique(pts_mean, 'rows');
    pts_mean_new = pts_mean(ia, :);
    tri = delaunay(pts_mean_new);
    MORPHED_IMAGE = morph(view_near, view_far, features(:, 1:2), features(:, 3:4), tri, 1-ratio, 1-ratio);
end

function [cost_i, cost_g, cost_f] = CalculateCost(near_idx, far_idx, im_near_subs, im_far_subs, windowsize)
im_near_sub = im_near_subs(:, :, near_idx);
im_far_sub = im_far_subs(:, :, far_idx);

intensity_near = im_near_sub(:);
intensity_far = im_far_sub(:);
intensity_cost = (intensity_near - intensity_far) .^ 2;
intensity_cost = sum(sum(intensity_cost));

[ nearGx, nearGy ] = imgradientxy(im_near_sub);
gradient_near = [ nearGx(:); nearGy(:) ];
[ farGx, farGy ] = imgradientxy(im_far_sub);
gradient_far = [ farGx(:); farGy(:) ];
gradient_cost = (gradient_near - gradient_far) .^ 2;
gradient_cost = sum(sum(gradient_cost));

cellSize = [windowsize/2 windowsize/2];
hog_near = extractHOGFeatures(im_near_sub, 'CellSize', cellSize);
hog_far = extractHOGFeatures(im_far_sub, 'CellSize', cellSize);
hog_cost = (hog_near - hog_far) .^ 2;
hog_cost = sum(sum(hog_cost));

cost_i = intensity_cost;
cost_g = gradient_cost;
cost_f = hog_cost;

% disp([num2str(cost_i) ', ' num2str(cost_g) ', ' num2str(cost_f)]);
% subplot(1, 2, 1);
% imshow(im_near_sub);
% subplot(1, 2, 2);
% imshow(im_far_sub);
% pause;
end

function IMAGE_MAT = inter8_mat(LF, P_r, P_1, P_2, U_r, U_1, U_2, H_r, H_1, H_2, c)
    global PARAMS;

    P_1(P_r == 1) = 0;    P_2(P_r == 0) = 0;
    U_1(U_r == 1) = 0;    U_2(U_r == 0) = 0;
    H_1(H_r == 1) = 0;    H_2(H_r == 0) = 0;
    
    if(c == 1)
        IMAGE_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 1)) + ...
                                       H_r .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 1))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 1)) + ...
                                       H_r .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 1)))))) + ...
                        ((P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 1)) + ...
                                       H_r .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 1))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 1)) + ...
                                       H_r .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 1))))));
    elseif(c == 2)
        IMAGE_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 2))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 2)))))) + ...
                        ((P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 2))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 2))))));
    elseif(c == 3)
        IMAGE_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 3))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 3)))))) + ...
                        ((P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_1 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 3))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_2) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + U_2 * (PARAMS.RESIZE_HEIGHT * 3) + H_2 * 3 + 3))))));
    end
end

function IMAGE = LoadFrame(LF, FRAME_NUM)
    global PARAMS;
    
    [TX, TY] = meshgrid(1:PARAMS.RESIZE_WIDTH, 1:PARAMS.RESIZE_HEIGHT);
    TX = TX(:);
    TY = TY(:);
    IMAGE = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
    IMAGE_TEMP = zeros(size(TX, 1), 3);

    IMAGE_TEMP(:, 1) = LF((FRAME_NUM-1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 3); 
    IMAGE_TEMP(:, 2) = LF((FRAME_NUM-1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 2);
    IMAGE_TEMP(:, 3) = LF((FRAME_NUM-1) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 1);

    IMAGE(:, :, 1) = reshape(IMAGE_TEMP(:, 1), [PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH]);
    IMAGE(:, :, 2) = reshape(IMAGE_TEMP(:, 2), [PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH]);
    IMAGE(:, :, 3) = reshape(IMAGE_TEMP(:, 3), [PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH]);
end

