clc
clear all;
close all;

A = imread('w1.png');
B = imread('2.png');

A=rgb2gray(A);
% A=imrotate(A,20);
imshow((A))
for i=1:size(A,1)
    for j=1:size(A,2)
        if A(i,j)>120
            A(i,j)=255;
        end
    end
end
imshow(A);
% %% quad 600um
% Rmin = 130;
% Rmax = 140;
% pitch=273;
%% hex 1000um
Rmin = 210;
Rmax = 240;
pitch=454;
[centersBright, radiiBright] = imfindcircles(A,[Rmin Rmax],'Method','TwoStage','Sensitivity',0.99);
% [centersDark, radiiDark] = imfindcircles(A,[Rmin Rmax],'ObjectPolarity','dark');
% viscircles(centersBright(1,:), radiiBright(1),'Color','r');
% viscircles(centersDark, radiiDark,'Color','r');

centr=(max(centersBright)-min(centersBright))/2+min(centersBright);

dist=centersBright-centr;
dist1=sqrt(dist(:,1).^2.+dist(:,2).^2);
[min,ind]=min(dist1);
viscircles(centersBright(ind,:), radiiBright(ind),'Color','r');
% viscircles(centersBright, radiiBright,'Color','r');





%% 
% refofusing algorithm
C=centersBright(ind,:);

pt=pitch/2;
pixels=20; %%num tot of pixeles to be overlapped
distances=centersBright-C;
normal_distances=(distances/pitch);
normal_distances1=round(normal_distances);
A=uint16(B);
I=uint16(zeros(size(A)));

% Video setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vid_name=('refoc_hex');
movie = VideoWriter(vid_name,'Uncompressed AVI');
movie.FrameRate = 5;
open(movie);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=0;
a=tic;
for pixel=-50-10:49-10
    count=count+1
    tic
    I=uint16(zeros(size(A)));
    for i=1:size(normal_distances,1)
        
% % %         I((C(1)-(pitch/2+pixels)):(C(1)+(pitch/2+pixels)), C(2)-(pitch/2+pixels):+(pitch/2+pixels))=I(C(1)-(pitch/2+pixels):C(1)+(pitch/2+pixels), C(2)-(pitch/2+pixels):+(pitch/2+pixels))...
% % %             +A(C(1)+normal_distances(1)*pitch+pixel*normal_distances(1)-(pitch/2+pixels):C(1)+normal_distances(1)*pitch+pixel*normal_distances(1)+(pitch/2+pixels), C(2)+normal_distances(2)*pitch+pixel*normal_distances(2)-(pitch/2+pixels):C(2)+normal_distances(2)*pitch+pixel*normal_distances(2)+(pitch/2+pixels));
% % % %                 I=I+A;
I(round(C(2)-(pitch/2+pixels)):round(C(2)+(pitch/2+pixels)), round(C(1)-(pitch/2+pixels)):round(C(1)+(pitch/2+pixels)),:)=I(round(C(2)-(pitch/2+pixels)):round(C(2)+(pitch/2+pixels)), round(C(1)-(pitch/2+pixels)):round(C(1)+(pitch/2+pixels)),:)...
            +A(round(C(2)+normal_distances(i,2)*pitch+pixel*normal_distances(i,2)-(pitch/2+pixels)):round(C(2)+normal_distances(i,2)*pitch+pixel*normal_distances(i,2)+(pitch/2+pixels)), ...
        round(C(1)+normal_distances(i,1)*pitch+pixel*normal_distances(i,1)-(pitch/2+pixels)):round(C(1)+normal_distances(i,1)*pitch+pixel*normal_distances(i,1)+(pitch/2+pixels)),:);
%                 I=I+A;
%     figure;
%     imshow(A(C(2)+normal_distances(i,2)*pitch+pixel*normal_distances(i,2)-(pitch/2+pixels):C(2)+normal_distances(i,2)*pitch+pixel*normal_distances(i,2)+(pitch/2+pixels), ...
%         C(1)+normal_distances(i,1)*pitch+pixel*normal_distances(i,1)-(pitch/2+pixels):C(1)+normal_distances(i,1)*pitch+pixel*normal_distances(i,1)+(pitch/2+pixels)));
% figure; imshow(uint8(I));

    end
    I1=I(round(C(2)-(pitch/2+pixels)):round(C(2)+(pitch/2+pixels)), round(C(1)-(pitch/2+pixels)):round(C(1)+(pitch/2+pixels)),:);
    I1=uint8(I1./size(normal_distances,1));
    
    frame = im2frame((I1));
    writeVideo(movie,frame);
%     figure; imshow(I);
toc
end
b=tic-a
close(movie);