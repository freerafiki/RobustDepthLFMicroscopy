clc
clear all;
close all;

filename = '/Users/Palma/Documents/Valencia/Images/Telescope/3.png';
filename1='../Images/Micro_1/FocalStack1/pitch473/1';
%read the image
I=double(imread((filename)));


pitch=473; % 682
pt=(pitch/2);
%center of EI(1)
C0=[295, 536];  %478, 831	
%EI in each row
a=[3 4 3];
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
%save(EI_naming);
%% refocusing part

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

%done for selecting the "center" of refocus
for k=5:5 %EI_number
    %refocusing video details
    text=[filename1 '_' num2str(k)];
    %movie = VideoWriter(text ,'Uncompressed AVI');
    %movie.FrameRate = 2;
    %open(movie);

    %pick up the several elemental images
    for step=ini:fin
        delta=step-pitch;
        offset=70; %%%focus offset of the optical system
        step_pix=13.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel
        posiz=-(step_pix*delta-offset);
        raccolta(1:EI_number)=struct('imag',Ii);
        for kk=1:EI_number
            C_new_im=int16([C_im(kk, 1)+(EI_naming(k,1)-EI_naming(kk,1))*step*cos(pi/6),C_im(kk, 2)+(EI_naming(k,2)-EI_naming(kk,2))/2*step]);
            raccolta(kk).imag((C_new_im(1)-pt):(C_new_im(1)+pt), (C_new_im(2)-pt):(C_new_im(2)+pt),:)= I((C_im(kk, 1)-pt):(C_im(kk, 1)+pt), (C_im(kk, 2)-pt):(C_im(kk, 2)+pt),:);
        end

    %here I am summing up the different images
    H=zeros(size(raccolta(1).imag));
    for n=1:EI_number
        H=(H+raccolta(n).imag);
    end
    %here I normalize the final image refocused at a certain plane
    H=uint8(H/EI_number);

    a=size(H);
    H1=H(round(C_im(k,1)-pt):round(C_im(k,1)+pt),round(C_im(k,2)-pt):round(C_im(k,2)+pt),:);

    a1=size(H1);
        for r=1:a1(1)
            for m=1:a1(2)
                if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-5)^2
                    H1(round(r), round(m),:)=0;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% STORE THE REFOCUSED IMAGES (H1) IN HERE %%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ref = strcat(filename1, '_center=', num2str(k), '_f=', num2str(posiz), '.bmp');       
        imwrite(H1, Ref);
        %frame = im2frame((H1));
        %writeVideo(movie,frame);
    end
    
        Icentral=raccolta(k).imag(round(C_im(k,1)-pt):round(C_im(k,1)+pt),round(C_im(k,2)-pt):round(C_im(k,2)+pt),:);
        a1=size(Icentral);
        for r=1:a1(1)
            for m=1:a1(2)
                if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-5)^2
                    Icentral(r,m,:)=0;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% STORE THE CENTRAL IMAGE (Icentral) IN HERE %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ref2 = strcat(filename1, '_center=', num2str(k), '.bmp');
        imwrite(uint8(Icentral), Ref2);
        
    %close(movie);
end