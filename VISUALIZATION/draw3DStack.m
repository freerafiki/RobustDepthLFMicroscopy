function h = draw3DStack(cv, matte)

if size(size(cv),2) == 4
    new_cv = zeros(size(cv,1), size(cv,2), size(cv, 4));
    for i = 1: size(cv,4)
        new_cv(:,:,i) = rgb2gray(uint8(cv(:,:,:,i))).*uint8(matte);
    end
else
    new_cv = cv;
end
[x,y,z] = meshgrid(1:size(new_cv,1), 1:size(new_cv,2), 1:size(new_cv,3));
xslice = [];
yslice = [];
zslice = linspace(1,size(new_cv,3),size(new_cv,3));
h = slice(x,y,z,new_cv,xslice,yslice,zslice);
hold on;
set(h,'LineStyle','none');