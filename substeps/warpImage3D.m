function I2warped = warpImage3D(I2, u, v, w)
%------------------------------------------------------------------------------
% Warp image based on given vectorfield by bilinear interpolation
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
I2=double(I2);
[height, width, length] = size(I2);

u1 = u;
u1(isnan(u1)) = 0;
v1 = v;
v1(isnan(v1)) = 0;
w1 = w;
w1(isnan(w1)) = 0;

[X, Y, Z]=meshgrid(1:width,1:height, 1:length);
 
xNew = X + u1;
yNew = Y + v1;
zNew = Z + w1;
xNew(xNew > width) = width;
xNew(xNew < 1) = 1;
yNew(yNew > height) = height;
yNew(yNew < 1) = 1;
zNew(zNew > length) = length;
zNew(zNew < 1) = 1;

I2warped=interp3(X, Y, Z, I2, xNew, yNew, zNew, 'bilinear');
I2warped(isnan(u)|isnan(v)|isnan(w)) = NaN;

end
