function L = lift_hv_nodetail(A000, A110, A011, A101, clift, which)
%------------------------------------------------------------------------------
% Apply hv lifting step with no detail stored 
%
% Warning: Currently only min lifting in this step is implemented
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------
tp = class(A000);

zD010 = zeros(size(A000, 1), size(A110, 2), size(A000, 3));
zD100 = zeros(size(A110, 1), size(A000, 2), size(A000, 3));
zD001 = zeros(size(A000, 1), size(A000, 2), size(A011, 3));
zD111 = zeros(size(A101, 1), size(A011, 2), size(A101, 3));

if strcmp(which,'min')

    A000 = A000 - min(zeros(size(A000),tp),...
        hvLift000min(zD010, zD100, zD001, clift));     % X1
    A110 = A110 - min(zeros(size(A110),tp), ...
        hvLift110min(zD100, zD010, zD111, clift));     % X1
    A011 = A011 - min(zeros(size(A011),tp), ...
        hvLift011min(zD001, zD111, zD010, clift));     % X1
    A101 = A101 - min(zeros(size(A101),tp), ...
        hvLift101min(zD111, zD001, zD100, clift));     % X1
   
%
%  Stage: reverse predict 1 
    A010 = zD010 + hvLift010min(A000, A110, A011, clift);                % Y1
    A100 = zD100 + hvLift100min(A110, A000, A101, clift);                % Y1
    A001 = zD001 + hvLift001min(A011, A000, A101, clift);                % Y1
    A111 = zD111 + hvLift111min(A101, A011, A110, clift);

   sizeR = size(A000) + size(A111);
%  Merge
    L = putcolor000(A000, sizeR) + putcolor101(A101, sizeR) + ...
        putcolor110(A110, sizeR) + putcolor011(A011, sizeR) + ...
        putcolor010(A010, sizeR) + putcolor100(A100, sizeR) + ...
        putcolor001(A001, sizeR) + putcolor111(A111, sizeR);
    
elseif strcmp(which,'max')
        error('not implemented')
        
else
    error('unknown lifting')
end