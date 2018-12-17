function [imgs, c_img] = read_imgs(path)

fprintf('\nReading the focal stack and the central image..\n');
I=double(imread((path)));
pitch=383; %473; % 682
pt=(pitch/2);
%center of EI(1)
C0=[344, 496]; %[295, 536];  %478, 831	
%EI in each row
a=[2 3 2];
%refocusing ini and fin
ini=445; %670;
fin=450; %694;
%number of rows
Nrow=size(a);
Nrow=Nrow(2);
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

Ii=zeros(size(I));
raccolta(1:EI_number)=struct('imag',Ii);
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

for k=1:EI_number
    %pick up the several elemental images
    for step=ini:fin
        delta=step-pitch;
        offset=70; %%%focus offset of the optical system
        step_pix=13.5 %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel
        posiz=-(step_pix*delta-offset);
        raccolta(1:EI_number)=struct('imag',Ii);
        for kk=1:EI_number
            C_new_im=int16([C_im(kk, 1)+(EI_naming(k,1)-EI_naming(kk,1))*step*cos(pi/6),C_im(kk, 2)+(EI_naming(k,2)-EI_naming(kk,2))/2*step]);
            raccolta(kk).imag((C_new_im(1)-pt):(C_new_im(1)+pt), (C_new_im(2)-pt):(C_new_im(2)+pt),:)= I((C_im(kk, 1)-pt):(C_im(kk, 1)+pt), (C_im(kk, 2)-pt):(C_im(kk, 2)+pt),:);
        end
    end
end