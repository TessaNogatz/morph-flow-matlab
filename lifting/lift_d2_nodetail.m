function [A000, A101] = lift_d2_nodetail(A000, yA101, zA101, clift, which)
%------------------------------------------------------------------------------
% Apply d2 lifting step with no detail stored 
%
% Warning: Currently only min lifting in this step is implemented
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------
tp = class(A000);
zD101 = zeros(yA101, size(A000, 2), zA101, tp);

zA000 = zeros(size(A000),tp);

if strcmp(which,'min')

    A000 = A000 -...
        min(zA000,d2Lift000min(zD101, size(A000, 1), size(A000, 3), clift));   
    A101 = zD101 + ...
              d2Lift101min(A000, size(zD101, 1), size(zD101, 3), clift);
          
elseif strcmp(which,'max')
    error('not implemented')
        
else
    error('unknown lifting')
end
