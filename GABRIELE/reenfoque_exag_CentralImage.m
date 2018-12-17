clc
clear all;
close all;
filename='/Users/Palma/Documents/Valencia/Images/Chip/1.png';
filename1='1_test';

%read the image
I=uint32(imread((filename)));

%defining the number of microlenses (Ncentr= in the center, Ntop= in the
%top
Ncentr=3;
Ntop=2;
pitch=486;

%overlapping of pixels that might have "focused" zones (ini=starting point,
%fin= end point
ini=497;
fin=510;

%scaling factor of the image (leave it to 1)
scale=1;

I=imresize(I,scale,'nearest');
ini1=(pitch-(pitch-ini))*scale;
fin1=(pitch-(pitch-fin))*scale;
pitch=pitch*scale;


%refocusing video details
movie = VideoWriter(filename1,'Uncompressed AVI');
movie.FrameRate = 2;
open(movie);
Ii=zeros(size(I));

%simply half of the pitch I'll need it later
pt=pitch/2;

%I am defining the number of microlenses array( ex. if I have an array of
%2-3-2 microlenses it means I have 3 horizontal arrays of microlenses
Liv=Ncentr-Ntop;
NumL=Liv*2+1;
num_images=0;

%this is delta Y I'm rounding it because it gives errors otherwise
Dy=round(pitch*cos(pi/6));
%raw and col are needed if the central image is not well centered in the
%whole image (it shoul actually be fine but just in case)
raw=0;
col=0;
%CN is the definition of the center of the central image
CN=[round(size(I,1)/2+raw),round((size(I,2)/2)+col)];
%counting how many images I expect (considering also black images or call
%them dummy images 
for n=Ntop:Ncentr
    if n==Ncentr
    num_images=num_images+n;
    else
    num_images=num_images+n*2;
    end
end


for step=ini1:fin1
    
    %%%%%%%%%%
    %It's calculating the depth of the focus at a certain number of pixel
    %overlapped
    delta=step-pitch;
    offset=70;
    step_pix=27.5/scale;
    posiz=-(step_pix*delta-offset);
    %%%%%%%%%%
    
    vistas=step;
    
    %raccolta is a struct that I will use to store all the images
    %separately
    raccolta(1:num_images)=struct('imag',Ii);
    
    h=-Liv;
    im_num=1;
    
    
    %% this is to make a struct where I define the center of the different images so that I can extract the sinlge view and store it in one of the component of "raccolta"
for n=Ntop:Ncentr
   
   if mod(n,2)==1 
       l=-(n-1)/2;
       for t=1:n
       C=[round(CN(1)+round(h*Dy)),round(CN(2)+l*pitch)];
       C1=[round(round(CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1)))), round(round(CN(2)+l*pitch+step/pitch*(CN(2)-C(2))))];
       raccolta(im_num).imag((C1(1)-pt):(C1(1)+pt), (C1(2)-pt):(C1(2)+pt),:)= I((C(1)-pt):(C(1)+pt), (C(2)-pt):(C(2)+pt),:);
       im_num=im_num+1;
       l=l+1;
       end
   else
       if mod(n,2)==0
           l=-n/2+0.5;
           for t=1:n
           C=[CN(1)+round(h*Dy), CN(2)+l*pitch];
           C1=[round(CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1))), round(CN(2)+l*pitch+step/pitch*(CN(2)-C(2)))];
           raccolta(im_num).imag(C1(1)-pt:C1(1)+pt, C1(2)-pt:C1(2)+pt,:)= I(C(1)-pt:C(1)+pt, C(2)-pt:C(2)+pt,:);
           im_num=im_num+1;
           l=l+1;
           end
       end
   end
   h=h+1;
end
n=0;
for n=Ncentr-1:-1:Ntop
   
   if mod(n,2)==1
       l=-(n-1)/2;
       for t=1:n
       C=[CN(1)+round(h*Dy), CN(2)+l*pitch];
       C1=[CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1)), CN(2)+l*pitch+step/pitch*(CN(2)-C(2))];
       raccolta(im_num).imag(round(C1(1)-pt):round(C1(1)+pt), round(C1(2)-pt):round(C1(2)+pt),:)= I(round(C(1)-pt):round(C(1)+pt), round(C(2)-pt):round(C(2)+pt),:);
       im_num=im_num+1;
       l=l+1;
       end
   else
       if mod(n,2)==0
           l=-n/2+0.5;
           for t=1:n
           C=[CN(1)+round(h*Dy), CN(2)+l*pitch];
           C1=[round(CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1))), round(CN(2)+l*pitch+step/pitch*(CN(2)-C(2)))];
           raccolta(im_num).imag(round(C1(1)-pt):round(C1(1)+pt), round(C1(2)-pt):round(C1(2)+pt),:)= I(round(C(1)-pt):round(C(1)+pt), round(C(2)-pt):round(C(2)+pt),:);
           im_num=im_num+1;
           l=l+1;
           end
       end
   end
   h=h+1;
end

%%
%here I am summing up the different images
H=zeros(size(raccolta(1).imag));
for n=1:num_images
    H=(H+raccolta(n).imag);
end
%here I normalize the final image refocused at a certain plane
H=uint8(H/num_images);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%HERE SHOULD GO THE DEPTH ESTIMATION CODE%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




a=size(H);
H1=H(round(a(1)/2-(pitch/2)):round(a(1)/2+(pitch/2)),round(a(2)/2-(pitch/2)):round(a(2)/2+(pitch/2)),:);

a1=size(H1);
    for r=1:a1(1)
        for m=1:a1(2)
            if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-5)^2
                H1(r,m,:)=0;
            end
        end
    end

text=['delta z=' num2str(posiz) '  '];
H1 = insertText(H1, [0 (a1(1)-30)], text, 'FontSize', 18, ...
           'BoxOpacity', 1,'BoxColor', 'black','TextColor', 'white');

imshow(H1);

end

%% get the central image
Icentral=raccolta(4).imag(round(a(1)/2-(pitch/2)):round(a(1)/2+(pitch/2)),round(a(2)/2-(pitch/2)):round(a(2)/2+(pitch/2)),:);
a1=size(Icentral);
    for r=1:a1(1)
        for m=1:a1(2)
            if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-5)^2
                Icentral(r,m,:)=0;
            end
        end
    end
Ref = strcat(filename1, '_CentralImage', '.bmp');
% Ref = strcat('_CentralImage', '.bmp');
imwrite(uint8(Icentral), Ref);       
%%
close(movie);