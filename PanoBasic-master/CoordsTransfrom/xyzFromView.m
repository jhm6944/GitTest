function [ out3DNorm, out3DPlane ] = xyzFromView( xy, xyz0, imgH, imgW )
%XYZFROMVIEW Project 2D image coordinates on a perspective image to 3D 
%points on panorama
%   INPUT:
%   xy: 2D image coordinates
%   xyz0: 3D position of the center of the perspective view. It equals to
%   v*R, where v is the viewing direction of the center of perspective
%   view, and R is the focal length.
%   imH,imW: the size of the perspective view.
%   OUTPUT:
%   out3DPlane: the 3D coordinates of those points on image plane
%   out3DNorm: normalized out3DPlane in unit vector

% im is the tangent plane, contacting with ball at [x0 y0 z0]
x0 = xyz0(1);
y0 = xyz0(2);
z0 = xyz0(3);

uv = xyz2uvN(xyz0,1);
x = uv(1);  y = uv(2);

vecposX = [cos(x) -sin(x) 0];
vecposY = cross([x0 y0 z0], vecposX);
vecposY = vecposY ./ norm(vecposY);

numPoint = size(xy,1);
out3DPlane = repmat([x0 y0 z0], numPoint, 1) ...
           + repmat(xy(:,1)-(imgW+1)/2, 1, 3).*repmat(vecposX, numPoint, 1) ...
           + repmat(xy(:,2)-(imgH+1)/2, 1, 3).*repmat(vecposY, numPoint, 1);

out3DNorm = out3DPlane ./ repmat(sqrt(sum(out3DPlane.^2, 2)), 1, 3);

end

