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

%     morph_image = renderingLF_morph(LF{1, 1}, LF{2, 1}, posx_front, posy_front);
    morph_image = renderingLF_morph_load(LF{1, 1}, LF{2, 1}, posx_front, posy_front);
    imshow(morph_image);
% end

