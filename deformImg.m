function [imgDef, uGT, vGT, wGT] = deformImg(imgOrig, maxDisp)
    [height, width, length] = size(imgOrig);
    [X, Y, Z] = meshgrid(1:width,1:height, 1:length);
    uGT = zeros(size(X));
    vGT = zeros(size(X));
    wGT = zeros(size(X));
    xNew = X;
    xNew(Z < round(0.5*length)) = X(Z < round(0.5*length)) + (Z(Z < round(0.5*length))-length)*0.005*maxDisp;
    uGT(Z < round(0.5*length)) = (Z(Z < round(0.5*length))-length)*0.005*maxDisp;
    xNew(Z >= round(0.5*length)) = X (Z >= round(0.5*length))+ Z(Z >= round(0.5*length))*0.005*maxDisp;
    uGT(Z >= round(0.5*length)) = Z(Z >= round(0.5*length))*0.005*maxDisp;
    yNew = Y;
    zNew = Z - 0.2 -(maxDisp-0.2)./sqrt(1+0.5.*exp(-0.05*(Z-round(length/2))));
    wGT = - 0.2 -(maxDisp-0.2)./sqrt(1+0.5.*exp(-0.05*(Z-round(length/2))));
    xNew(xNew > width) = width;
    xNew(xNew < 1) = 1;
    yNew(yNew > height) = height;
    yNew(yNew < 1) = 1;
    zNew(zNew > length) = length;
    zNew(zNew < 1) = 1;

    imgDef=interp3(X, Y, Z, double(imgOrig), xNew, yNew, zNew, 'linear');
    
    