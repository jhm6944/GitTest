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
        
        features_temp = cat(1, features_temp, [w, h_n, w, h, cost(y, x)]);
    end
    [C,ia,ic] = unique(features_temp(:, 1:2), 'rows');
    features = cat(1, features, features_temp(ia, :));
end

view_near = padarray(view_near, [25, 25], 0, 'both');
view_far = padarray(view_far, [25, 25], 0, 'both');

for i=1:size(features, 1)
    w_n = features(i, 1);
    h_n = features(i, 2);
    w_f = features(i, 3);
    h_f = features(i, 4);
    
    view_near_sub = view_near((h_n-20):(h_n+20), (w_n-20):(w_n+20), :);
    view_far_sub = view_far((h_f-20):(h_f+20), (w_f-20):(w_f+20), :);
    
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





