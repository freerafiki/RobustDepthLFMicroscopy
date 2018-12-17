function [imgs, c_img] = read_EI(path)

fprintf('\nReading the focal stack and the central image..\n');
I=double(imread((path)));
pitch=683;%telescope444; %473; % 682
pt=floor(pitch/2);
C0=[478, 831];%[295, 536];%telescope[344, 496];
a=[2 3 2];
Nrow = size(a,2);
EI_number=sum(a);
EI_naming=zeros(EI_number,2);
b=1;
if a(1)<a(2)
    f=2;
else
    f=1;
end
for k=1:Nrow
    l=0;
    for t=1:a(k)
        EI_naming(b,:)=[k,f+l];
        l=l+2;
        b=b+1;
    end
    if f==1
        f=2;
    else
        f=1;        
    end
end

C_im=zeros(size(EI_naming));
%defining centers of the original EIs
%pick up the several elemental images
if a(1)<a(2)
        for k1=1:EI_number
            C_im(k1,:)=[C0(1)+(cos(pi/6)*pitch*(EI_naming(k1,1)-1)), C0(2)+pitch/2*(EI_naming(k1,2)-2)];
        end
else
        for k1=1:EI_number
            C_im(k1,:)=[C0(1)+(cos(pi/6)*pitch*(EI_naming(k1,1)-1)), C0(2)+pitch/2*(EI_naming(k1,2)-1)];
        end
end

for kk=1:EI_number
    %pick up the several elemental images
    C_new_im=int16([C_im(kk, 1)+(EI_naming(kk,1))*cos(pi/6),C_im(kk, 2)+(EI_naming(kk,2))/2]);
    imgs{kk} = I((C_new_im(1)-pt):(C_new_im(1)+pt), (C_new_im(2)-pt):(C_new_im(2)+pt),:);
    %imshow(uint8(imgs{kk}));
end
c_k = floor(EI_number/2);
c_img = imgs{c_k};