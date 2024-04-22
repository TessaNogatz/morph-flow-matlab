function [dx, dy, dz] = gradient110hv(A110, A000, A011, A101, cmin)
%------------------------------------------------------------------------------
% Calculate gradient for 110 in dodecahedral setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 
[m110, n110, l110] = size(A110);
[m000, n000, l000] = size(A000);
[m101, n101, l101] = size(A101);
[m011, n011, l011] = size(A011);

 
cellVol = 16 / 9 * sqrt(2) ^ 2 * sqrt(3);
adjustedFaceArea = 2 / 3 ; 
% To the actual face area formula * 2 * sqrt(2) is missing. This is
% cancelled out with the factors from unit normal computation

% rhomboid and filled

if n110 == n000 % n even
    if m110 == m000 % m even
        sUeD = stripU3D(extD3D(A000, cmin));
        sLeR = stripL3D(extR3D(A000, cmin));
        R1 = - A000 + stripL3D(extR3D(sUeD, cmin));
        R2x = + sLeR - sUeD;
        R2y = - R2x;
    elseif m110 == m000-1 % m odd
        sD = stripD3D(A000);
        sU = stripU3D(A000);
        sDsLeR = stripL3D(extR3D(sD, cmin));
        sLeRsU = stripL3D(extR3D(sU, cmin));
        R1 = - sD + sLeRsU;
        R2x = sDsLeR - sU;
        R2y = -R2x;
    else
        error(' gradient110hv - sizes do not match ');
    end
elseif n110 == n000-1 % n odd
    if m110 == m000 % m even
        sR = stripR3D(A000);
        sL = stripL3D(A000);
        sUeD = stripU3D(extD3D(A000, cmin));
        R1 = -sR + stripL3D(sUeD);
        R2x = sL - stripR3D(sUeD);
        R2y = -R2x;
    elseif m110 == m000-1 % m odd
        sD = stripD3D(A000);
        sU = stripU3D(A000);
        sRsD = stripR3D(sD);
        sLsU = stripL3D(sU);
        sLsD = stripL3D(sD);
        sRsU = stripR3D(sU);
        R1 = sRsD + sLsU;
        R2x = + sLsD - sRsU;
        R2y = - R2x;
    else
        error(' gradient110hv - sizes do not match ');
    end
else 
    error(' gradient110hv - sizes do not match ');
end

% rhomboid and triangle

if m110 == m011 % n even
    if l110 == l011 % l even
        sUeD = stripU3D(extD3D(A011, cmin));
        sBeF = stripB3D(extF3D(A011, cmin));
        sBeFsUeD = stripB3D(extF3D(sUeD, cmin));
        S1y = -A011 + sBeFsUeD;
        S1z = -S1y;
        S2 = -sBeF + sUeD;
    elseif l110 == l011+1 % l odd
        eB = extB3D(A011, cmin);
        eF = extF3D(A011, cmin);
        sUeD = stripU3D(extD3D(A011, cmin));
        eFsUeD = extF3D(sUeD, cmin);
        eBsUeD = extB3D(sUeD, cmin);
        S1y = - eB + eFsUeD;
        S1z = - S1y;
        S2 = -eF + eBsUeD;
    else
        error(' gradient110hv - sizes do not match ');
    end
elseif m110 == m011-1 % n odd
    if l110 == l011 % m even
        sD = stripD3D(A011); 
        sU = stripU3D(A011);
        sBeF = stripB3D(extF3D(A011, cmin));
        S1y = -sD + stripU3D(sBeF);
        S1z = -S1y;
        S2 = - stripD3D(sBeF) + sU;
    elseif l110 == l011+1 % m odd
        sU = stripU3D(A011);
        sD = stripD3D(A011);
        S1y = -extB3D(sD, cmin) + extF3D(sU, cmin);
        S1z = -S1y;
        S2 = -extF3D(sD, cmin) + extB3D(sU, cmin);
    else
        error(' gradient110hv - sizes do not match ');
    end
else 
    error(' gradient110hv - sizes do not match ');
end

% rhomboid and star

if n110 == n101 % n even
    if l110 == l101 % l even
        sBeF = stripB3D(extF3D(A101, cmin));
        sLeR = stripL3D(extR3D(A101, cmin));
        sLeRsBeF = stripL3D(extR3D(sBeF, cmin));
        T1x = -A101 + sLeRsBeF;
        T1z = -T1x;
        T2 = sLeR - sBeF;
    elseif l110 == l101+1 % l odd
        eB = extB3D(A101, cmin);
        eF = extF3D(A101, cmin);
        sLeR = stripL3D(extR3D(A101, cmin));
        T1x = -eB + extF3D(sLeR, cmin);
        T1z = -T1x;
        T2 = extB3D(sLeR, cmin) - eF;
    else
        error(' gradient110hv - sizes do not match ');
    end
elseif n110 == n101-1 % n odd
    if l110 == l101 % m even
        sR = stripR3D(A101);
        sL = stripL3D(A101);
        sBeF = stripB3D(extF3D(A101, cmin));
        T1x = -sR + stripL3D(sBeF);
        T1z = -T1x;
        T2 = + sL - stripR3D(sBeF);
    elseif l110 == l101+1 % m odd
        eB = extB3D(A101, cmin);
        eF = extF3D(A101, cmin);
        sReB = stripR3D(eB);
        sLeF = stripL3D(eF);
        sLeB = stripL3D(eB);
        sReF = stripR3D(eF);
        T1x = -sReB + sLeF;
        T1z = -T1x;
        T2 = sLeB + sReF;
    else
        error(' gradient110hv - sizes do not match ');
    end
else 
    error(' gradient110hv - sizes do not match ');
end

dx = adjustedFaceArea / cellVol * ( R1 + R2x + T1x + T2);
dy = adjustedFaceArea / cellVol * ( R1 + R2y + S1y + S2);
dz = adjustedFaceArea / cellVol * ( S1z + S2 + T1z + T2);
