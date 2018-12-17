function trimap_out = labelExpansion(img_input,cdata)
% Expand trimap by geodesic distance

I = double(img_input);
trimap = cdata;
grayI = double(rgb2gray(img_input));
J = stdfilt(grayI,ones(7,7));
J = ordfilt2(J,1,ones(9,9));

% computing geodesic distance
M = double(trimap < 255);
DF = mexGetGeodesicDis(I,J,M);

M = double(trimap > 0);
DB = mexGetGeodesicDis(I,J,M);

% computing threshold
MF = (cdata == 255);
MB = (cdata == 0);
for k=1:1:100
    SE = [0,1,0;1,1,1;0,1,0];
    MF = imdilate(MF,SE);
    if max(MF(:) & MB(:)) == 1
        break;
    end
end
threshold = min((k-1)/2,20);
%threshold = 20;

% get new trimap
EF = DF < threshold;   EB = DB < threshold;
trimap(EF>0) = 255;     trimap(EB>0) = 0;
trimap_out = trimap;

end
