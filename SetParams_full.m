function LF = SetParams_full()

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
LF_FILE{1, 1} = ['Row7'; 'Row8'; 'Row9'];
LF_FILE{2, 1} = ['Row13'; 'Row14'; 'Row15'];

LF_FILE{3, 1} = ['Column166'; 'Column167'; 'Column168'];
LF_FILE{4, 1} = ['Column110'; 'Column111'; 'Column112'];

% LF_FILE{1, 1} = ['Row15'; 'Row16'; 'Row17'];
% LF_FILE{2, 1} = ['Row21'; 'Row22'; 'Row23'];

for i=1:2
    LF{i, 1} = zeros(3*PARAMS.LFU_W * PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3, 1, 'uint8');
    for n=1:PARAMS.LFU_W*3
        if(n <= PARAMS.LFU_W)
            file = sprintf('E:/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(1, :), n);
        elseif(n <= PARAMS.LFU_W*2)
            file = sprintf('E:/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(2, :), n - PARAMS.LFU_W);
        else
            file = sprintf('E:/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(3, :), n - PARAMS.LFU_W*2);
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

for i=3:4
    LF{i, 1} = zeros(3*PARAMS.LFU_W * PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3, 1, 'uint8');
    for n=1:PARAMS.LFU_W*3
        if(n <= PARAMS.LFU_W)
            file = sprintf('E:/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(1, :), n);
        elseif(n <= PARAMS.LFU_W*2)
            file = sprintf('E:/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(2, :), n - PARAMS.LFU_W);
        else
            file = sprintf('E:/BMW_5K/BMW_ALIGN/%s/%04d.jpg' , LF_FILE{i, 1}(3, :), n - PARAMS.LFU_W*2);
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
        nn = (PARAMS.LFU_W * 3 - n);
        
        LF{i, 1}((nn) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF{i, 1}((nn) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF{i, 1}((nn) * (PARAMS.RESIZE_HEIGHT * PARAMS.RESIZE_WIDTH * 3) + (TX - 1) * (PARAMS.RESIZE_HEIGHT * 3) + (TY - 1) * 3 + 1) = img_C3(:);
        
        disp(['LF ' num2str(i) ' - ' num2str(n) ', ' num2str(nn)]);
    end
end

end

