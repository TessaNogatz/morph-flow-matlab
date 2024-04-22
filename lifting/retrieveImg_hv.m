function A000 = retrieveImg_hv(C, S, lev, A000, A110, A011, A101, clift, which)
%------------------------------------------------------------------------------
% Utility function to perform image retrival from hv setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------  
tp = class(A000);

if strcmp(which,'min')
    [D010, D100, D001, D111] = retrievehvLift(lev-2, 'd', C, S);
       A000 = A000 - min(zeros(size(A000),tp),...
            hvLift000min(D010, D100, D001, clift));     % X1
       A110 = A110 - min(zeros(size(A110),tp), ...
            hvLift110min(D100, D010, D111, clift));     % X1
       A011 = A011 - min(zeros(size(A011),tp), ...
            hvLift011min(D001, D111, D010, clift));     % X1
       A101 = A101 - min(zeros(size(A101),tp), ...
            hvLift101min(D111, D001, D100, clift));     % X1

    %
    %  Stage: reverse predict 1 
       A010 = D010 + hvLift010min(A000, A110, A011, clift);                % Y1
       A100 = D100 + hvLift100min(A110, A000, A101, clift);                % Y1
       A001 = D001 + hvLift001min(A011, A000, A101, clift);                % Y1
       A111 = D111 + hvLift111min(A101, A011, A110, clift);
       clear D010 D100 D001 D111;

       %sizeR = size(A000) + size(A110) + size(A011) + size(A101) + ...
       %       size(A010) + size(A100) + size(A001) + size(A111);
       sizeR = size(A000) + size(A111);
    %  Merge
       A000 = putcolor000(A000, sizeR) + putcolor101(A101, sizeR) + ...
              putcolor110(A110, sizeR) + putcolor011(A011, sizeR) + ...
              putcolor010(A010, sizeR) + putcolor100(A100, sizeR) + ...
              putcolor001(A001, sizeR) + putcolor111(A111, sizeR);
          
 elseif strcmp(which,'max')
    error('not implemented')
        
else
    error('unknown lifting')
end
