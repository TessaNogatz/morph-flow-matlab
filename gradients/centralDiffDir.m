function [Ix, Iy, Iz] = centralDiffDir(I)
%------------------------------------------------------------------------------
% Calculate central gradients with 5-stencil report
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------
k2 = (1/12)*[1 -8 0 8 -1];

Ix = imfilter(I, k2, 'replicate');
Iy = imfilter(I, reshape(k2, 5, 1, 1), 'replicate');
Iz = imfilter(I, reshape(k2, 1, 1, 5), 'replicate');



end