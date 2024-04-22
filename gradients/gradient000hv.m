function [dx, dy, dz] = gradient000hv(A000, A110, A011, A101, cmin)
%------------------------------------------------------------------------------
% Calculate gradient for 000 in dodecahedral setting
%
% Design and implementation in 3D
% (c) 2024 Dr. Tessa Nogatz, tessa.nogatz@gmail.com
%------------------------------------------------------------------------------ 
[m000, n000, l000] = size(A000);
[m110, n110, l110] = size(A110);
[m011, n011, l011] = size(A011);
[m101, n101, l101] = size(A101);
 
cellVol = 16 / 9 * sqrt(2) ^ 2 * sqrt(3);
adjustedFaceArea = 2 / 3 ; 
% To the actual face area formula * 2 * sqrt(2) is missing. This is
% cancelled out with the factors from unit normal computation

% filled and rhombiod
if n000 == n110 % n even
    if m000 == m110 % m even
        sDeU = stripD3D(extU3D(A110, cmin));
        eLsRsDeU = extL3D(stripR3D(sDeU), cmin);
        eLsR = extL3D(stripR3D(A110), cmin);
        S1  = -eLsRsDeU + A110;
        S2x = - eLsR + sDeU;
        S2y = -S2x;
    elseif m000 == m110+1 % m odd
        eU = extU3D(A110, cmin);
        eD = extD3D(A110, cmin);
        eLsReD = extL3D(stripR3D(eD), cmin);
        S1  = -extL3D(stripR3D(eU), cmin) + eD;
        S2x = - eLsReD + eU;
        S2y = -S2x;
    else
        error(' gradient000hv - sizes do not match ');
    end
elseif n000 == n110+1 % n odd
    if m000 == m110 % m even
        eR = extR3D(A110, cmin);
        eL = extL3D(A110, cmin);
        sDeUeR = stripD3D(extU3D(eR, cmin));
        S1  = - extU3D(stripD3D(eL), cmin) + eR;
        S2x = - eL + sDeUeR;
        S2y = -S2x;
    elseif m000 == m110+1 % m odd
        eL = extL3D(A110, cmin);
        eR = extR3D(A110, cmin);
        eLeD = extD3D(eL, cmin);
        eUeR = extU3D(eR, cmin);
        S1  = - extU3D(eL, cmin) + extD3D(eR, cmin);
        S2x = - eLeD + eUeR;
        S2y = -S2x;
    else
        error(' gradient000hv - sizes do not match ');
    end
else 
    error(' gradient000hv - sizes do not match ');
end

% filled and triangle
if n000 == n011 % n even
    if l000 == l011 % l even
        eLsR = extL3D(stripR3D(A011), cmin);
        sBeF = stripB3D(extF3D(A011, cmin));
        R1 = extF3D(stripB3D(eLsR), cmin) + A011;
        R2x = -eLsR + sBeF;
        R2z = -R2x;
    elseif l000 == l011+1 % l odd
        eB = extB3D(A011, cmin);
        eF = extF3D(A011, cmin);
        eBeLsR = extL3D(stripR3D(eB), cmin);
        R1 = - extL3D(stripR3D(eF), cmin) + eB;
        R2x = - eBeLsR + eF;
        R2z = -R2x;
    else
        error(' gradient000hv - sizes do not match ');
    end
elseif n000 == n011+1 % n odd
    if l000 == l011 % m even
        eL = extL3D(A011, cmin);
        eR = extR3D(A011, cmin);
        eRsBeF = stripB3D(extF3D(eR, cmin));
        R1 = extF3D(stripB3D(eL), cmin) + eR;
        R2x = - eL + eRsBeF;
        R2z = -R2x;
    elseif l000 == l011+1 % m odd
       eF = extF3D(A011, cmin);
       eB = extB3D(A011, cmin);
       eLeB = extL3D(eB, cmin);
       eReF = extR3D(eF, cmin);
       R1 = -extL3D(eF, cmin) + extR3D(eB, cmin);
       R2x = - eLeB + eReF;
       R2z = -R2x;
    else
        error(' gradient000hv - sizes do not match ');
    end
else 
    error(' gradient000hv - sizes do not match ');
end

% filled and star

if l000 == l101 % n even
        if m000 == m101 % m even
            eFsB = extF3D(stripB3D(A101), cmin);
            sDeU = stripD3D(extU3D(A101, cmin));
            T2z = - eFsB + sDeU;
            T2y = -T2z;
            T1 = + A101 - extU3D(stripD3D(eFsB), cmin);
        elseif m000 == m101+1 % m odd
            eU = extU3D(A101, cmin);
            eD = extD3D(A101, cmin);
            eFsBeD = extF3D(stripB3D(eD), cmin);
            T2z = - eFsBeD + eU;
            T2y = -T2z;
            T1  = eD - extF3D(stripB3D(eU), cmin);
        else
        error(' gradient000hv - sizes do not match ');
        end
elseif l000 == l101+1 % n odd
        if m000 == m101 % m even
            eB = extB3D(A101, cmin);
            eF = extF3D(A101, cmin);
            eUsD = extU3D(stripD3D(A101), cmin);
            eFeUsD = extF3D(eUsD, cmin);
            T2z = - eB + eFeUsD;
            T2y = -T2z;
            T1 = eF - extB3D(eUsD, cmin);
        elseif m000 == m101+1 % m odd
            eD = extD3D(A101, cmin);
            eU = extU3D(A101, cmin);
            eFeD = extF3D(eD, cmin);
            eBeU = extB3D(eU, cmin);
            T2z = - eFeD + eBeU;
            T2y = -T2z;
            T1 = extB3D(eD, cmin) + extF3D(eU, cmin);
        else
            error(' gradient000hv - sizes do not match ');
        end
else 
    error(' gradient000hv - sizes do not match ');
end

dx = adjustedFaceArea / cellVol * ( S1 + S2x + R1 + R2x);
dy = adjustedFaceArea / cellVol * ( S1 + S2y + T1 + T2y);
dz = adjustedFaceArea / cellVol * ( R1 + R2z + T1 + T2z); % R1x - R2 - T1 - T2z

