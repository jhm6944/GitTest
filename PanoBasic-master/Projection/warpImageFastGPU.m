function im_warp = warpImageFastGPU(im, XXdense, YYdense)

%{
Citation:
J. Xiao, K. A. Ehinger, A. Oliva and A. Torralba.
Recognizing Scene Viewpoint using Panoramic Place Representation.
Proceedings of 25th IEEE Conference on Computer Vision and Pattern Recognition, 2012.
http://sun360.mit.edu
%}

minX = max(1,floor(min(min(XXdense)))-1);
minY = max(1,floor(min(min(YYdense)))-1);

maxX = min(size(im,2),ceil(max(max(XXdense)))+1);
maxY = min(size(im,1),ceil(max(max(YYdense)))+1);

% im = im(minY:maxY,minX:maxX,:);
im1 = gpuArray(im(minY:maxY,minX:maxX,1));
im2 = gpuArray(im(minY:maxY,minX:maxX,2));
im3 = gpuArray(im(minY:maxY,minX:maxX,3));

im_warp(:,:,1) = interp2(im1, XXdense-minX+1, YYdense-minY+1,'LINEAR');
im_warp(:,:,2) = interp2(im2, XXdense-minX+1, YYdense-minY+1,'LINEAR');
im_warp(:,:,3) = interp2(im3, XXdense-minX+1, YYdense-minY+1,'LINEAR');
% 
% for c=1:size(im,3)
%     % im_warp(:,:,c) = uint8(interp2(double(im(:,:,c)), XXdense-minX+1, YYdense-minY+1,'*cubic'));
%     %im_warp(:,:,c) = interp2(im(:,:,c), XXdense-minX+1, YYdense-minY+1,'*cubic');
%     im_warp(:,:,c) = interp2(im(:,:,c), XXdense-minX+1, YYdense-minY+1,'*linear');
% %     im_warp(:,:,c) = interp2(im(:,:,c), XXdense-minX+1, YYdense-minY+1,'*nearest');
% %     im_warp(:,:,c) = interp2(im(:,:,c), XXdense-minX+1, YYdense-minY+1,'linear');
% end
