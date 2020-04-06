%% Set Params
addpath('PanoBasic-master');
add_path;

global PARAMS;

PARAMS.LFU_W = 50;
PARAMS.HEIGHT = 2560;
PARAMS.WIDTH  = 5120;
PARAMS.TIMES = 500;

PARAMS.HEIGHT_RESIZE = 512;
PARAMS.WIDTH_RESIZE  = 1024;

PARAMS.IMG_SIZE = 31;
PARAMS.IMG_FOV = pi()/12;
PARAMS.IMG_STD_TH = 0.04;
PARAMS.IMG_COST_TH = 50;

[TX, TY] = meshgrid(1:PARAMS.WIDTH_RESIZE, 1:PARAMS.HEIGHT_RESIZE);
TX = TX(:);
TY = TY(:);

LF = cell(2, 1);
LF_FILE = cell(2, 1);
LF_FILE{1, 1} = ['Row7' ; 'Row8' ; 'Row9' ];
LF_FILE{2, 1} = ['Row13'; 'Row14'; 'Row15'];

for i=1:length(LF)
    LF{i, 1} = zeros(3*PARAMS.LFU_W * PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3, 1, 'uint8');
    for n=1:PARAMS.LFU_W*3
        if(n <= PARAMS.LFU_W)
            file = sprintf('../Matlab_FeatureExtraction/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(1, :), n);
        elseif(n <= PARAMS.LFU_W*2)
            file = sprintf('../Matlab_FeatureExtraction/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(2, :), n - PARAMS.LFU_W);
        else
            file = sprintf('../Matlab_FeatureExtraction/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(3, :), n - PARAMS.LFU_W*2);
        end
        img = imread(file);
        img = imresize(img, [512, 1024]);
        img_C1 = img(:, :, 1);
        img_C2 = img(:, :, 2);
        img_C3 = img(:, :, 3);
        
%         if(i == 2)
%             nn = (PARAMS.LFU_W * 3 - n);
%         else
%             nn = n - 1;
%         end
        nn = n - 1;
        
        LF{i, 1}((nn) * (PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3) + (TX - 1) * (PARAMS.HEIGHT_RESIZE * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF{i, 1}((nn) * (PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3) + (TX - 1) * (PARAMS.HEIGHT_RESIZE * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF{i, 1}((nn) * (PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3) + (TX - 1) * (PARAMS.HEIGHT_RESIZE * 3) + (TY - 1) * 3 + 1) = img_C3(:);
        
        disp(['LF ' num2str(i) ' - ' num2str(n) ', ' num2str(nn)]);
    end
end

%% DEFINE PARAMETERS
ANGLE_S = -1 * deg2rad(45);
ANGLE_E = deg2rad(45);
OUT_W = floor((ANGLE_E - ANGLE_S) / deg2rad(360/PARAMS.WIDTH_RESIZE));
W=(1:OUT_W)';
ANGLE = ANGLE_S + deg2rad(360/PARAMS.WIDTH_RESIZE) .* (W-1);

U = ANGLE.*(180.0 / pi)*(1.0 / 180.0)*PARAMS.WIDTH_RESIZE / 2 + PARAMS.WIDTH_RESIZE / 2;
U = round(U);
U(U < 1) = U(U < 1) + PARAMS.WIDTH_RESIZE;
U(U > PARAMS.WIDTH_RESIZE) = U(U > PARAMS.WIDTH_RESIZE) - PARAMS.WIDTH_RESIZE;

%% FIND FEATURES
height = PARAMS.HEIGHT_RESIZE;
width = PARAMS.WIDTH_RESIZE;

[TX, TY] = meshgrid(1:width, 1:height);
LON = (TX/(width/2) - 1) .* pi();
LAT = -1*(TY/(height/2) - 1) .* pi()/2;

windowsize = 16;
fov = pi()/12;

disparity = zeros(PARAMS.LFU_W*3, 2*height, width);
cost = zeros(PARAMS.LFU_W*3, 2*height, width);
for n=1:PARAMS.LFU_W*3
    im_near = rgb2gray(LoadFrame(LF{1, 1}, n)); 
    P = -1 * (PARAMS.LFU_W) .* tan(ANGLE) + n;% - PARAMS.LFU_W/2 + n + (PARAMS.LFU_W * 3)/2;
    P = round(P);
    
    for idx_w = 1:size(P, 1)
        if(P(idx_w, 1) >= 1 && P(idx_w, 1) <= (PARAMS.LFU_W*3))
            im_far = rgb2gray(LoadFrame(LF{2, 1}, P(idx_w, 1)));
            
            im_near_subs = zeros(windowsize, windowsize, 2*height);
            im_far_subs = zeros(windowsize, windowsize, 2*height);
            disparityLine = zeros(2*height, 1);
            costLine = zeros(2*height, 1);
            
            % Preload sub-images
            tic;
            w = U(idx_w, 1);
            parfor idx_h=1:2*height
                h = idx_h - height/2;
                w_n = 0;
                if(h < 1)
                    h = -1 * h;
                    w_n = w + width/2;
                elseif(h > height)
                    h = 2 * height - h;
                    w_n = w + width/2;
                end
                
                if(w_n > width)
                    w_n = w_n - width;
                end
                
                lon = (w_n/(width/2) - 1) * pi();
                lat = -1*(h/(height/2) - 1) * pi()/2;
                
%                 disp([num2str(w) ', ' num2str(h) ', ' num2str(lon) ', ' num2str(lat)]);
                
                im_near_sub = imgLookAt(im2double(im_near), lon, lat, windowsize, fov);
                im_near_sub(isnan(im_near_sub)) = 0;
                im_near_subs(:, :, idx_h) = im_near_sub;
                
                im_far_sub = imgLookAt(im2double(im_far), lon, lat, windowsize, fov);
                im_far_sub(isnan(im_far_sub)) = 0;
                im_far_subs(:, :, idx_h) = im_far_sub;
                
%                 subplot(1, 2, 1);
%                 imshow(im_near_sub);
%                 subplot(1, 2, 2);
%                 imshow(im_far_sub);
%                 pause;
            end
            toc;
            
            % Calculate Costs
            bound = 40;
            tic;
%             fp = fopen('cost_debug.txt', 'w');
            parfor farY=1:2*height
                cost_buff = zeros(bound, 3);
                count = 0;
                if(farY < height)
                    for nearY = farY:-1:max(farY-bound, 1)
                        count = count + 1;  
                        [cost_i, cost_g, cost_f] = CalculateCost(nearY, farY, im_near_subs, im_far_subs, windowsize);
                        cost_buff(count, :) = [cost_i cost_g cost_f];
%                         fprintf(fp, '%f %f %f\n', cost_i, cost_g, cost_f);
                    end
                    cost_t = cost_buff(1:count, :);
                    cost_t2 = cost_t(:, 1) + cost_t(:, 2) / 15 + 50 * cost_t(:, 3);
                    [val, loc] = min(cost_t2);
                    disparityLine(farY, 1) = loc;
                    costLine(farY, 1) = val;
                elseif(farY > height)
                    for nearY = farY:min(farY+bound, 2*height)
                        count = count + 1;
                        [cost_i, cost_g, cost_f] = CalculateCost(nearY, farY, im_near_subs, im_far_subs, windowsize);
                        cost_buff(count, :) = [cost_i cost_g cost_f];
%                         fprintf(fp, '%f %f %f\n', cost_i, cost_g, cost_f);
                    end
                    cost_t = cost_buff(1:count, :);
                    cost_t2 = 0.25 * cost_t(:, 1) + 0.25 * cost_t(:, 2) + 0.5 * cost_t(:, 3);
                    [val, loc] = min(cost_t2);
                    disparityLine(farY, 1) = loc;
                    costLine(farY, 1) = val;
                else
                    disparityLine(farY, 1) = 0;
                    costLine(farY, 1) = 0;
                end
            end
            toc;
%             fclose(fp);
            disparity(n, :, idx_w) = disparityLine;
            cost(n, :, idx_w) = costLine;
        else
            disparity(n, :, idx_w) = -1;
            cost(n, :, idx_w) = 9999;
        end
    end
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
end

function IMAGE = LoadFrame(LF, FRAME_NUM)
global PARAMS;

[TX, TY] = meshgrid(1:PARAMS.WIDTH_RESIZE, 1:PARAMS.HEIGHT_RESIZE);
TX = TX(:);
TY = TY(:);
IMAGE = zeros(PARAMS.HEIGHT_RESIZE, PARAMS.WIDTH_RESIZE, 3, 'uint8');
IMAGE_TEMP = zeros(size(TX, 1), 3);

IMAGE_TEMP(:, 1) = LF((FRAME_NUM-1) * (PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3) + (TX - 1) * (PARAMS.HEIGHT_RESIZE * 3) + (TY - 1) * 3 + 3); 
IMAGE_TEMP(:, 2) = LF((FRAME_NUM-1) * (PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3) + (TX - 1) * (PARAMS.HEIGHT_RESIZE * 3) + (TY - 1) * 3 + 2);
IMAGE_TEMP(:, 3) = LF((FRAME_NUM-1) * (PARAMS.HEIGHT_RESIZE * PARAMS.WIDTH_RESIZE * 3) + (TX - 1) * (PARAMS.HEIGHT_RESIZE * 3) + (TY - 1) * 3 + 1);

IMAGE(:, :, 1) = reshape(IMAGE_TEMP(:, 1), [PARAMS.HEIGHT_RESIZE, PARAMS.WIDTH_RESIZE]);
IMAGE(:, :, 2) = reshape(IMAGE_TEMP(:, 2), [PARAMS.HEIGHT_RESIZE, PARAMS.WIDTH_RESIZE]);
IMAGE(:, :, 3) = reshape(IMAGE_TEMP(:, 3), [PARAMS.HEIGHT_RESIZE, PARAMS.WIDTH_RESIZE]);
end
