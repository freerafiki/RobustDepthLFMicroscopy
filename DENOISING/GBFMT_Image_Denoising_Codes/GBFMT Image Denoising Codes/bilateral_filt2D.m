function [yb yg]=bilateral_filt2D(x,sigmas,sigmar,ksize)

%%% 2D Bilateral Filter.

%%% Author : B.K. SHREYAMSHA KUMAR
%%% Created on 06-05-2011.
%%% Updated on 06-05-2011. 

half_ksize=floor(ksize/2);

%%% Gaussian Kernel Generation. 
gk=gauss_ker2D(sigmas,ksize);

%%% To take care of boundaries.
[M,N]=size(x);
xm=zeros(M+ksize-1,N+ksize-1);
xm(half_ksize+1:M+half_ksize,half_ksize+1:N+half_ksize)=x;
xm(1:half_ksize,:)=xm(ksize:-1:half_ksize+2,:);	%% Row-wise periodic extension.
xm(M+half_ksize+1:M+ksize-1,:)=xm(M+half_ksize-1:-1:M,:);
xm(:,1:half_ksize)=xm(:,ksize:-1:half_ksize+2);	%% Column-wise periodic extension.
xm(:,N+half_ksize+1:N+ksize-1)=xm(:,N+half_ksize-1:-1:N);

%%% Bilateral Filter Implementation.
[MM,NN]=size(xm);
for ii=half_ksize+1:MM-half_ksize
   for jj=half_ksize+1:NN-half_ksize
      xtemp=xm(ii-half_ksize:ii+half_ksize,jj-half_ksize:jj+half_ksize);
      pix_diff=abs(xtemp-xtemp(half_ksize+1,half_ksize+1));
      rgk=exp(-(pix_diff.^2/(2*(sigmar)^2)));
      yg(ii-half_ksize,jj-half_ksize)=sum(sum(xtemp.*gk))/sum(sum(gk));
      yb(ii-half_ksize,jj-half_ksize)=sum(sum(xtemp.*gk.*rgk))/sum(sum(gk.*rgk));
   end
end