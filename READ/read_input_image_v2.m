function [fs, map, EIs] = read_input_image(filename, type, pitch, ini, fin, a, C0, offset, step_pix)

%% GENERATE FOCAL STACK AND ELEMENTAL IMAGES
% input: path is the path to the image (FiMic image)
%
% output: 
%   - fs: matrices [h,w,n] with [h,w] the size of an image and n the number
%   of focal planes generated (images stack one on top of the other
%   - map: corresponding relations with depth and disparity. 
%   The m-th image (fs(:,:,m)) is refocused at dists(m,1) [actual depth] 
%   and is shifted of m pixels [actual disparity]).
%   - EIs: the elemental images 



%% GABRIELE's METHOD
fprintf('\nReading the input image and creating focal stack and elemental images..\n');
tic
%read image
I=double(imread((filename)));
%defining the number of microlenses (Ncentr= in the center, Ntop= in the
%top
Ncentr=3;
Ntop=2;
%pitch=455;


%overlapping of pixels that might have "focused" zones (ini=starting point,
%fin= end point
%ini=497;
%fin=510;


%% OUTPUTs
num_of_EIs = Ncentr+Ntop*2;
num_of_planes = fin - ini+1;
EIs = zeros(pitch, pitch, 3, num_of_EIs);
fs = zeros(pitch, pitch, 3, num_of_planes);
map = zeros(num_of_planes, 1);

%scaling factor of the image (leave it to 1)
scale=1;

I=imresize(I,scale,'nearest');
ini1=(pitch-(pitch-ini))*scale;
fin1=(pitch-(pitch-fin))*scale;
pitch=pitch*scale;


Ii=zeros(size(I));

%simply half of the pitch I'll need it later
pt=floor(pitch/2);

%I am defining the number of microlenses array( ex. if I have an array of
%2-3-2 microlenses it means I have 3 horizontal arrays of microlenses
Liv=Ncentr-Ntop;
NumL=Liv*2+1;
num_images=Ntop+Ncentr;

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
           C1=int16([round(round(CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1)))), round(round(CN(2)+l*pitch+step/pitch*(CN(2)-C(2))))]);
           raccolta(im_num).imag((C1(1)-pt):(C1(1)+pt), (C1(2)-pt):(C1(2)+pt),:)= I((C(1)-pt):(C(1)+pt), (C(2)-pt):(C(2)+pt),:);
           im_num=im_num+1;
           l=l+1;
           end
       else
           if mod(n,2)==0
               l=-n/2+0.5;
               for t=1:n
                   C=[CN(1)+round(h*Dy), CN(2)+l*pitch];
                   C1=int16([round(CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1))), round(round(CN(2)+l*pitch+step/pitch*(CN(2)-C(2))))]);
                   raccolta(im_num).imag((C1(1)-pt):(C1(1)+pt), (C1(2)-pt):(C1(2)+pt),:)= I((C(1)-pt):(C(1)+pt), (C(2)-pt):(C(2)+pt),:);
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
           C1=int16([CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1)), CN(2)+l*pitch+step/pitch*(CN(2)-C(2))]);
           raccolta(im_num).imag(round(C1(1)-pt):round(C1(1)+pt), round(C1(2)-pt):round(C1(2)+pt),:)= I(round(C(1)-pt):round(C(1)+pt), round(C(2)-pt):round(C(2)+pt),:);
           im_num=im_num+1;
           l=l+1;
           end
        else
           if mod(n,2)==0
               l=-n/2+0.5;
               for t=1:n
               C=[CN(1)+round(h*Dy), CN(2)+l*pitch];
               C1=int16([round(CN(1)+round(h*Dy)+step/pitch*(CN(1)-C(1))), round(CN(2)+l*pitch+step/pitch*(CN(2)-C(2)))]);
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
    H=uint8(H/7);
    
    a=size(H);
    H1=H(round(a(1)/2-floor(pitch/2)):round(a(1)/2+floor(pitch/2)),round(a(2)/2-floor(pitch/2)):round(a(2)/2+floor(pitch/2)),:);

    a1=size(H1);
    for r=1:a1(1)
        for m=1:a1(2)
            if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-5)^2
                H1(r,m,:)=0;
            end
        end
    end

    %text=['delta z=' num2str(posiz) '  '];
    %H1 = insertText(H1, [0 (a1(1)-30)], text, 'FontSize', 18, ...
    %           'BoxOpacity', 1,'BoxColor', 'black','TextColor', 'white');

    current_index = step - ini1 + 1;
    fprintf('Focal Stack creation: focal plane %d of %d\n', current_index, fin-ini+1);
    if size(size(H1),2) == 3
        fs(:,:,:,current_index) = H1;        
    elseif size(size(H1),2) == 2
        fs(:,:,1,current_index) = H1;
        fs(:,:,2,current_index) = H1;
        fs(:,:,3,current_index) = H1;
    else
        fprintf('Wrong dimensions of refocused image - something went wrong');
    end
    map(current_index) = posiz; % actual depth of the focal plane
end

%% get the central image
for c = 1:num_of_EIs
    Icentral=raccolta(c).imag(round(a(1)/2-floor(pitch/2)):round(a(1)/2+floor(pitch/2)),round(a(2)/2-floor(pitch/2)):round(a(2)/2+floor(pitch/2)),:);
    a1=size(Icentral);
    for r=1:a1(1)
        for m=1:a1(2)
            if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-5)^2
                Icentral(r,m,:)=0;
            end
        end
    end
    EIs(:,:,:,c) = Icentral;
end
    
fprintf('Done!                       ');
toc
