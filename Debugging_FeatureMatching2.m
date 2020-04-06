LF_A0 = LF{1, 1};
LF_A1 = LF{2, 1};
POS_X = 0;
POS_Y = 0;

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

view_near = cat(1, rot90(VIEW_BACK_NEAR(1:size(VIEW_BACK_NEAR, 1)/2, :, :), 2), ...
    VIEW_FRONT_NEAR, rot90(VIEW_BACK_NEAR((size(VIEW_BACK_NEAR, 1)/2 + 1):size(VIEW_BACK_NEAR, 1), :, :), 2));

% view_near = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
% view_near(:, (PARAMS.RESIZE_WIDTH/2 - PARAMS.RESIZE_WIDTH/8 + 1):(PARAMS.RESIZE_WIDTH/2 + PARAMS.RESIZE_WIDTH/8), :) = VIEW_FRONT_NEAR;
% view_near(:, 1:PARAMS.RESIZE_WIDTH/8, :) = VIEW_BACK_NEAR(:, (size(VIEW_BACK_NEAR, 2)/2 + 1):size(VIEW_BACK_NEAR, 2), :);
% view_near(:, (PARAMS.RESIZE_WIDTH - PARAMS.RESIZE_WIDTH/8 + 1):PARAMS.RESIZE_WIDTH, :) = VIEW_BACK_NEAR(:, 1:size(VIEW_BACK_NEAR, 2)/2, :);

VIEW_FRONT_FAR_3 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 1));
VIEW_FRONT_FAR_2 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 2));
VIEW_FRONT_FAR_1 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_FRONT-1, 0, 0, H-1, 0, 3));
VIEW_FRONT_FAR = cat(3, VIEW_FRONT_FAR_1, VIEW_FRONT_FAR_2, VIEW_FRONT_FAR_3);

VIEW_BACK_FAR_3 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 1));
VIEW_BACK_FAR_2 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 2));
VIEW_BACK_FAR_1 = im2uint8(inter8_mat(LF_A1, 0, P_far-1, 0, 0, U_BACK-1, 0, 0, H-1, 0, 3));
VIEW_BACK_FAR = cat(3, VIEW_BACK_FAR_1, VIEW_BACK_FAR_2, VIEW_BACK_FAR_3);

view_far = cat(1, rot90(VIEW_BACK_FAR(1:size(VIEW_BACK_FAR, 1)/2, :, :), 2), ...
    VIEW_FRONT_FAR, rot90(VIEW_BACK_FAR((size(VIEW_BACK_FAR, 1)/2 + 1):size(VIEW_BACK_FAR, 1), :, :), 2));

% view_far = zeros(PARAMS.RESIZE_HEIGHT, PARAMS.RESIZE_WIDTH, 3, 'uint8');
% view_far(:, (PARAMS.RESIZE_WIDTH/2 - PARAMS.RESIZE_WIDTH/8 + 1):(PARAMS.RESIZE_WIDTH/2 + PARAMS.RESIZE_WIDTH/8), :) = VIEW_FRONT_FAR;
% view_far(:, 1:PARAMS.RESIZE_WIDTH/8, :) = VIEW_BACK_FAR(:, (size(VIEW_BACK_FAR, 2)/2 + 1):size(VIEW_BACK_FAR, 2), :);
% view_far(:, (PARAMS.RESIZE_WIDTH - PARAMS.RESIZE_WIDTH/8 + 1):PARAMS.RESIZE_WIDTH, :) = VIEW_BACK_FAR(:, 1:size(VIEW_BACK_FAR, 2)/2, :);

height = PARAMS.HEIGHT/4;
width = PARAMS.WIDTH/4;
h_t = ((height*2)-1) / floor(height*2/10);
hh = 1:h_t:(height*2);
w_t = (OUT_W-1) / floor(OUT_W/10);
ww = 1:w_t:OUT_W;

load('disparity');
load('cost');
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

view_near = padarray(view_near, [25, 25], 0, 'both');
view_far = padarray(view_far, [25, 25], 0, 'both');

for i=1:size(features, 1)
    w_n = round(features(i, 1)) + 25;
    h_n = round(features(i, 2)) + 25;
    w_f = round(features(i, 3)) + 25;
    h_f = round(features(i, 4)) + 25;
    
    view_near_sub = view_near((h_n-25):(h_n+25), (w_n-25):(w_n+25), :);
    view_far_sub = view_far((h_f-25):(h_f+25), (w_f-25):(w_f+25), :);
    
    disp([num2str(w_n) ', ' num2str(h_n) ', ' num2str(w_f) ', ' num2str(h_f) ]);
    subplot(1, 2, 1);
    imshow(view_near_sub);
    subplot(1, 2, 2);
    imshow(view_far_sub);
    pause;
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





