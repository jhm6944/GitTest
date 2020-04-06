%% DEBUG FEATURES
global PARAMS;

load('disparity');
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