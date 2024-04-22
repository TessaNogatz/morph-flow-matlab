function [A000, A110, A011, A101] = lift_d1_nodetail(A000, A101, xD110, clift, which)
%------------------------------------------------------------------------------
% Apply d1 lifting step with no detail stored 
%
% Warning: Currently only max lifting in this step is implemented
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 
tp = class(A000);
zD110 = zeros(size(A101, 1), xD110, size(A000, 3), tp);
zD011 = zeros(size(A000, 1), xD110, size(A101, 3), tp);

zVF000 = zeros(size(A000),tp);
zVF101 = zeros(size(A101),tp);
para.sigma = 1;

if strcmp(which,'max')
    A000 = A000 - ...
          max(zVF000, d1Lift000max(zD110, zD011, size(A000, 2), clift));
    A101 = A101 - ...
          max(zVF101, d1Lift101max(zD110, zD011, size(A101, 2), clift));

    A110 = zD110 + d1Lift110max(A000, A101, size(zD110, 2), clift);
    A011 = zD011 + d1Lift011max(A000, A101, size(zD011, 2) ,clift);
elseif strcmp(which,'min')
        error('not implemented')
        
else
    error('unknown lifting')
end