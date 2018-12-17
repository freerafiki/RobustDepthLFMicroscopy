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
%%
fprintf('\nReading the input image and creating focal stack and elemental images..\n');
tic
%read image
I=double(imread((filename)));

pt=floor(pitch/2);
%number of rows
Nrow=size(a);
Nrow=Nrow(2);
EI_number=sum(a);
EI_naming=zeros(EI_number,2);
num_of_planes = fin - ini+1;
num_of_EIs = sum(a);
%% OUTPUTs
EIs = zeros(pitch, pitch, 3, num_of_EIs);
fs = zeros(pitch, pitch, 3, num_of_planes);
map = zeros(num_of_planes, 1);


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
c_img_index = ceil(num_of_EIs/2);
raccolta(1:EI_number)=struct('imag',Ii);
%% FOCAL STACK BUILT ONLY FOR CENTRAL IMAGE
for k=1:EI_number %EI_number
    fprintf('EI number %d of %d..\n', k, EI_number);
    %pick up the several elemental images
    
    for step=ini:fin
        delta=step-pitch;

        posiz=-(step_pix*delta-offset);
        
        for kk=1:EI_number
            x_ = round(C_im(kk, 1)+(EI_naming(k,1)-EI_naming(kk,1))*step*cos(pi/6));
            y_ = round(C_im(kk, 2)+(EI_naming(k,2)-EI_naming(kk,2))*step*.5); %.5 = sen(pi/6);
            C_new_im=int16([x_,y_]);
            %fprintf('c_new_im %d,%d, pt %d, C_im(kk) %d, %d\n', C_new_im(1), C_new_im(2), pt, C_im(kk,1), C_im(kk,2));
            raccolta(kk).imag(round(C_new_im(1)-pt):round(C_new_im(1)+pt), round(C_new_im(2)-pt):round(C_new_im(2)+pt),:)= I(round(C_im(kk, 1)-pt):round(C_im(kk, 1)+pt), round(C_im(kk, 2)-pt):round(C_im(kk, 2)+pt),:);
        end
        if k==c_img_index

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
            current_index = step - ini + 1;
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
            % current index is also the disparity - the amount of pixels
            % of the shift to create this image
            map(current_index) = posiz; % actual depth of the focal plane
            %frame = im2frame((H1));
            %writeVideo(movie,frame);
        end
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
    if size(size(Icentral),2) == 3
        EIs(:,:,:,k) = Icentral;        
    elseif size(size(Icentral),2) == 2
        EIs(:,:,1,k) = Icentral; 
        EIs(:,:,2,k) = Icentral; 
        EIs(:,:,3,k) = Icentral; 
    else
        fprintf('Wrong dimensions of central image - something went wrong');
    end

end

fprintf('Done!                       ');
toc
