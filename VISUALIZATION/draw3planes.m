function [xy, xz, yz] = draw3planes(x,y,cv)

xy = cv(x,y,:);
[nothing, z] = min(xy);
xz = cv(x,:,:);
yz = cv(:,y,:);

subplot(131)
imagesc(xy)
subplot(132)
imagesc(xz)
subplot(133)
imagesc(yz)