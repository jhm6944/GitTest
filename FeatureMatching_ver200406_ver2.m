LF = SetParams();
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

%% Functions
function MORPHED_IMAGE = renderingLF_morph(LF_A0, LF_A1, POS_X, POS_Y)%, MatchedFeatures)
    global PARAMS;
    ANGLE_S = -1 * deg2rad(45);
    ANGLE_E = deg2rad(45);

    [VIEW_FRONT_NEAR, VIEW_FRONT_FAR, VIEW_BACK_NEAR, VIEW_BACK_FAR] = GenImages(LF_A0, LF_A1, POS_X, POS_Y);
    [view_near, view_far] = ConversionFormOfImage(VIEW_FRONT_NEAR, VIEW_FRONT_FAR, VIEW_BACK_NEAR, VIEW_BACK_FAR, 1);

    height = PARAMS.HEIGHT/4;
    width = PARAMS.WIDTH/4;
    
    windowsize = 50;
    fov_1 = pi()/12;

    OUT_W = floor((ANGLE_E - ANGLE_S) / deg2rad(360/width));
    h_t = ((height*2)-1) / floor(height*2/20);
    hh = 1:h_t:(height*2);
    w_t = (OUT_W-1) / floor(OUT_W/20);
    ww = 1:w_t:OUT_W;
    
    disparity = zeros(size(hh, 2), size(ww, 2));
    cost = zeros(size(hh, 2), size(ww, 2));
    
%     for idx_w=1:size(ww, 2)
%         disp(['idx_w : ' num2str(idx_w) '/' num2str(size(ww, 2))]);
%         tic; 
%         w = ww(1, idx_w);
%         im_near = padarray(rgb2gray(view_near), [windowsize/2, windowsize/2], 0, 'both');
%         im_far = padarray(rgb2gray(view_far), [windowsize/2, windowsize/2], 0, 'both');
%         
%         im_far_subs  = zeros(windowsize+1, windowsize+1, size(hh, 2), 'uint8');
%         for idx_h=1:size(hh, 2)
%             h = hh(1, idx_h);
%             w_n = round(w) + 25;
%             h_n = round(h) + 25;
% 
%            im_far_sub = im_far((h_n-25):(h_n+25), (w_n-25):(w_n+25));
%            im_far_subs(:, :, idx_h) = im_far_sub;
% %            imshow(im_far_sub);
% %            pause;
%         end
% 
%         %% Extract Near images
%         im_near_subs = zeros(windowsize+1, windowsize+1, 2*height, 'uint8');
%         for idx_h=1:2*height
%             h = idx_h;
%             w_n = round(w) + 25;
%             h_n = round(h) + 25;
%             
%            im_near_sub = im_near((h_n-25):(h_n+25), (w_n-25):(w_n+25));
%            im_near_subs(:, :, idx_h) = im_near_sub;
%         end
%        
%         
%         %% Calculate cost
%         disparityLine = zeros(size(hh, 2), 1);
%         costLine = zeros(size(hh, 2), 1);
%         bound = 50;
%         parfor idx_h_far=1:size(hh, 2)
%             farY = round(hh(1, idx_h_far));
%             
%             if(idx_h_far >= 90 && idx_h_far <= 104)
%                 disparityLine(idx_h_far, 1) = 1;
%                 costLine(idx_h_far, 1) = -1;
%             elseif(farY < height)
%                 cost_buff = zeros(bound, 3);
%                 count = 0;
%                 for idx_h_near=farY:-1:max(farY - bound, 1)
%                     count = count + 1;
%                     [cost_i, cost_g, cost_f] = CalculateCost(idx_h_near, idx_h_far, im_near_subs, im_far_subs, windowsize);
%                     cost_buff(count, :) = [cost_i cost_g cost_f];
%                 end
%                 cost_t = cost_buff(1:count, :);
% %                 cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 250 * cost_t(:, 3);
%                 cost_t2 = cost_t(:, 1) + 250 * cost_t(:, 3);                
%                 [val, loc] = min(cost_t2);
%                 disparityLine(idx_h_far, 1) = loc;
%                 costLine(idx_h_far, 1) = val;
%             elseif(farY > height)
%                 cost_buff = zeros(bound, 3);
%                 count = 0;
%                 for idx_h_near = farY:min(farY + bound, 2*height)
%                     count = count + 1;
%                     [cost_i, cost_g, cost_f] = CalculateCost(idx_h_near, idx_h_far, im_near_subs, im_far_subs, windowsize);
%                     cost_buff(count, :) = [cost_i cost_g cost_f];
%                 end
%                 cost_t = cost_buff(1:count, :);
% %                 cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 250 * cost_t(:, 3);
%                 cost_t2 = cost_t(:, 1) + 250 * cost_t(:, 3);
%                 [val, loc] = min(cost_t2);
%                 disparityLine(idx_h_far, 1) = loc;
%                 costLine(idx_h_far, 1) = val;
%             else
%                 disparityLine(idx_h_far, 1) = 1;
%                 costLine(idx_h_far, 1) = -1;
%             end
%         end
%          toc;
%         %%
%         
% %         for idx_h_far=1:size(hh, 2)
% %             farY = round(hh(1, idx_h_far));
% %             
% %             c = costLine(idx_h_far, 1);
% %             loc = disparityLine(idx_h_far, 1) - 1;
% %             disp(['h_idx : ' num2str(idx_h_far) ', ' num2str(loc) ', ' num2str(c)]);
% %             disp(['std  : ' num2str(std(double(im_near_sub(:)))) ', ' num2str(std(double(im_far_sub(:))))]);
% %             
% %             if(farY < height)
% %                 im_near_sub = im_near_subs(:, :, farY - loc);
% %                 im_far_sub = im_far_subs(:, :, idx_h_far);
% %             else
% %                 im_near_sub = im_near_subs(:, :, farY + loc);
% %                 im_far_sub = im_far_subs(:, :, idx_h_far);
% %             end
% %             subplot(1, 2, 1);
% %             imshow(im_near_sub);
% %             subplot(1, 2, 2);
% %             imshow(im_far_sub);
% %             pause;
% %         end
% 
%         disparity(:, idx_w) = disparityLine;
%         cost(:, idx_w) = costLine;
%     end
    
    %% Extract features
    load('disparity');
    load('cost');
    features = zeros(0, 5);
    for x=1:size(disparity, 2)    
        features_temp = zeros(0, 5);
        for y=1:size(disparity, 1)
            w = ww(1, x);
            h = hh(1, y);
            
            if(h < height)
                h_n = h - (disparity(y, x) - 1);
                %             h_n = hh(1, y - (disparity(y, x) - 1));
            elseif(h > height)
                h_n = h + (disparity(y, x) - 1);
                %             h_n = hh(1, y + (disparity(y, x) - 1));
            else
                h_n = h;
            end
            
            h_near = h_n - height/2;
            h_far = h - height/2;
            
            if(h_near < 1)
                w_near = OUT_W - w +1;
            elseif(h_near > height)
                w_near = OUT_W - w +1;
            end
            
            if(h_far < 1)
                w_far = OUT_W - w +1;
            elseif(h_far > height)
                w_far = OUT_W - w +1;
            end
            
            features_temp = cat(1, features_temp, [w_near, h_n, w_far, h, cost(y, x)]);
        end
        
        y=2;
        while (y ~= size(features_temp, 1) - 1)
            curr_y = features_temp(y, 2);
            if(y < height)
                prev_y_idx = y - 1;
            elseif(y > height)
                prev_y_idx = y + 1;
            end
            prev_y = features_temp(prev_y_idx, 2);
            if(curr_y < prev_y)
                features_temp(y, :) = [];
            elseif(curr_y == prev_y)
                if(features_temp(y, 5) > features_temp(prev_y_idx, 5))
                    features_temp(y, :) = [];
                else
                    features_temp(prev_y_idx, :) = [];
                end
            else
                y = y + 1;
            end
        end
        features = cat(1, features, features_temp);
    end
    
    [view_near, view_far] = ConversionFormOfImage(VIEW_FRONT_NEAR, VIEW_FRONT_FAR, VIEW_BACK_NEAR, VIEW_BACK_FAR, 1);
%     
%     for i=1:size(features, 1)
%         showMatchedFeatures(view_near, view_far, features(i, 1:2), features(i, 3:4));
%         pause;
%     end
    
%     showMatchedFeatures(view_near, view_far, features(:, 1:2), features(:, 3:4));
    %% Mophing images
    ratio = (POS_Y + 50)/100;
    pts_mean = (features(:, 1:2) + features(:, 3:4)) / 2;
    tri = delaunay(pts_mean);
    
    
    MORPHED_IMAGE = morph(view_near, view_far, features(:, 1:2), features(:, 3:4), tri, 1-ratio, 1-ratio);
end

function [VIEW_FRONT_NEAR, VIEW_FRONT_FAR, VIEW_BACK_NEAR, VIEW_BACK_FAR] = GenImages(LF_A0, LF_A1, POS_X, POS_Y)
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
end

function [view_near, view_far] = ConversionFormOfImage(VIEW_FRONT_NEAR, VIEW_FRONT_FAR, VIEW_BACK_NEAR, VIEW_BACK_FAR, mode)
global PARAMS;
if(mode == 0)
    view_near = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
    view_near(:, (PARAMS.RESIZE_WIDTH/2 - PARAMS.RESIZE_WIDTH/8 + 1):(PARAMS.RESIZE_WIDTH/2 + PARAMS.RESIZE_WIDTH/8), :) = VIEW_FRONT_NEAR;
    view_near(:, 1:PARAMS.RESIZE_WIDTH/8, :) = VIEW_BACK_NEAR(:, (size(VIEW_BACK_NEAR, 2)/2 + 1):size(VIEW_BACK_NEAR, 2), :);
    view_near(:, (PARAMS.RESIZE_WIDTH - PARAMS.RESIZE_WIDTH/8 + 1):PARAMS.RESIZE_WIDTH, :) = VIEW_BACK_NEAR(:, 1:size(VIEW_BACK_NEAR, 2)/2, :);
    
    view_far = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
    view_far(:, (PARAMS.RESIZE_WIDTH/2 - PARAMS.RESIZE_WIDTH/8 + 1):(PARAMS.RESIZE_WIDTH/2 + PARAMS.RESIZE_WIDTH/8), :) = VIEW_FRONT_FAR;
    view_far(:, 1:PARAMS.RESIZE_WIDTH/8, :) = VIEW_BACK_FAR(:, (size(VIEW_BACK_FAR, 2)/2 + 1):size(VIEW_BACK_FAR, 2), :);
    view_far(:, (PARAMS.RESIZE_WIDTH - PARAMS.RESIZE_WIDTH/8 + 1):PARAMS.RESIZE_WIDTH, :) = VIEW_BACK_FAR(:, 1:size(VIEW_BACK_FAR, 2)/2, :);
elseif(mode == 1)
    view_near = cat(1, rot90(VIEW_BACK_NEAR(1:size(VIEW_BACK_NEAR, 1)/2, :, :), 2), ...
        VIEW_FRONT_NEAR, rot90(VIEW_BACK_NEAR((size(VIEW_BACK_NEAR, 1)/2 + 1):size(VIEW_BACK_NEAR, 1), :, :), 2));
    view_far = cat(1, rot90(VIEW_BACK_FAR(1:size(VIEW_BACK_FAR, 1)/2, :, :), 2), ...
        VIEW_FRONT_FAR, rot90(VIEW_BACK_FAR((size(VIEW_BACK_FAR, 1)/2 + 1):size(VIEW_BACK_FAR, 1), :, :), 2));
end
end

function [cost_i, cost_g, cost_f] = CalculateCost(near_idx, far_idx, im_near_subs, im_far_subs, windowsize)
im_near_sub = im_near_subs(:, :, near_idx);
im_far_sub = im_far_subs(:, :, far_idx);

intensity_near = im2double(im_near_sub(:));
intensity_far = im2double(im_far_sub(:));
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

% disp(['cost : ' num2str(cost_i) ', ' num2str(cost_f) ', ' num2str(cost_i + 250 * cost_f)]);
% disp(['std  : ' num2str(std(double(im_near_sub(:)))) ', ' num2str(std(double(im_far_sub(:))))]);
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