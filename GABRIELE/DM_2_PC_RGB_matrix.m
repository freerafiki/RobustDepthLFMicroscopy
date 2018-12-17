clc
clear all;
close all;
cd '/Users/Palma/Documents/Valencia/Results/Chip19_2405_phi06/depths/'
pitch=455;
%EI in each row
a=[2 3 2];
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

filename=('img1_depth.png');
%read the images
I=double(imread((filename)));
I_stack=(zeros([size(I),4*EI_number]));
for i = 4 : 4 % which image is using
    filename=strcat('img',num2str(i+1),'_no_norm_depth','.png');
    filename_im=strcat('img',num2str(i+1),'.png');
    I_stack(:,:,(i*4+1):(i*4+3))=(imread((filename_im)));
    I_stack(:,:,4*(i+1))=double(imread((filename)));
%     I_stack(:,:,4*(i+1))=I_stack(:,:,4*(i+1))./(max(max(I_stack(:,:,4*(i+1))))).*29;

% metto zero raggio esterno
a1=size(I_stack);

    for r=1:a1(1)
        for m=1:a1(2)
            if (r-a1(1)/2)^2+(m-a1(2)/2)^2>=(pitch/2-20)^2
                    I_stack(r,m,(i*4+1):(i*4+4))=0;
            end
        end
    end


end

%%

cd '/Users/Palma/Documents/Valencia/Code/matlab_depth/GABRIELE/'
mex Write_PC_txt.cpp
cd ..
Write_PC_txt(I_stack,EI_naming);
