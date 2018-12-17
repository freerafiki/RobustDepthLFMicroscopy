 function FM = focusmeasure(Image, Measure, WSize)
%This function measures the relative degree of focus of 
%an image. It may be invoked as:
%
%   FM = focusmeasure(Image, Method, WSize)
%
%Where 
%   Image,  is a DOUBLE Image and FM is a 
%           matrix the same size as Image with the 
%           computed focus measure for every pixel. 
%   WSize,  is the size of the neighborhood used to 
%           compute the focus value of every pixel. 
%           If WSize = 0, a single focus measure is 
%           computed for the whole image and, in this 
%           case, FM is a scalar. 
%   Method, is the focus measure algorithm as a string.
%           
% S. Pertuz
% Jan/2016


AF = (WSize==0)||isempty(WSize)||~exist('WSize','var');
if AF, WSize = 15;
end
MEANF = fspecial('average',[WSize WSize]);
switch upper(Measure)
    case 'ACMO' % Absolute Central Moment (Shirvaikar2004)
        if AF, FM = AcMomentum(im2uint8(Image));
        else
        FM = nlfilter((Image), [WSize WSize], @AcMomentum);
        FM = imfilter(FM, MEANF,'replicate');        
        end
        
    case 'BREN' % Brenner's (Santos97)
        [M, N] = size(Image);
        DH = zeros(M, N);
        DV = zeros(M, N);
        DV(1:M-2,:) = Image(3:end,:)-Image(1:end-2,:);
        DH(:,1:N-2) = Image(:,3:end)-Image(:,1:end-2);
        FM = max(DH, DV);        
        FM = FM.^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF,'replicate');
        end
        
    case 'CHEB' % Chebyshev moments-based (Yap2004)
        %Image = im2double(Image);
        if AF, FM = TchebiFocus(Image);
        else
        FM = nlfilter(Image, [WSize WSize], @TchebiFocus);
        FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'CONT' % Image contrast (Nanda2001)
        ImContrast = @(x) sum(abs(x(:)-x(5)));
        FM = nlfilter(Image, [3 3], ImContrast);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
    case 'CURV' % Image Curvature (Helmli2001)       
        M1 = [-1 0 1;-1 0 1;-1 0 1];
        M2 = [1 0 1;1 0 1;1 0 1];
        P0 = imfilter(Image, M1, 'replicate', 'conv')/6;
        P1 = imfilter(Image, M1', 'replicate', 'conv')/6;
        P2 = 3*imfilter(Image, M2, 'replicate', 'conv')/10 ...
            -imfilter(Image, M2', 'replicate', 'conv')/5;
        P3 = -imfilter(Image, M2, 'replicate', 'conv')/5 ...
            +3*imfilter(Image, M2, 'replicate', 'conv')/10;
        FM = abs(P0) + abs(P1) + abs(P2) + abs(P3);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'DCTE' % DCT energy ratio (Shen2006)
        FM = nlfilter(Image, [8 8], @DctRatio);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'DCTR' % DCT reduced energy ratio (Lee2009)
        FM = nlfilter(Image, [8 8], @ReRatio);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'DCTM' % DCT Modified (Lee2008)
        M = [1 1 -1 -1;1 1 -1 -1;...
            -1 -1 1 1;-1 -1 1 1];
        FM = imfilter(Image, M, 'replicate');
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'EIGV' % Eigenvalues-based (Wee2007)
        %Image = im2double(Image);
        if AF, FM = EigenFocus(Image);
        else        
        FM = nlfilter(Image, [WSize WSize], @EigenFocus);            
        FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'GDER' % Gaussian derivative (Geusebroek2000)
        N = floor(WSize/2);
        sig = N/3;
        [x,y] = meshgrid(-N:N, -N:N);
        G = exp(-(x.^2+y.^2)/(2*sig^2))/(2*pi*sig);
        Gx = -x.*G/(sig^2);Gx = Gx/sum(Gx(:));
        Gy = -y.*G/(sig^2);Gy = Gy/sum(Gy(:));
        Rx = imfilter(Image, Gx, 'conv', 'replicate');
        Ry = imfilter(Image, Gy, 'conv', 'replicate');
        FM = Rx.^2+Ry.^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'GLVA' % Graylevel variance (Krotkov86)        
        if AF, FM = std2(Image).^2;
        else FM = stdfilt(Image, ones(WSize,WSize)).^2;
            FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'GLLV' %Graylevel local variance (Pech2000)
        if AF
            LVar = stdfilt(Image, ones(WSize,WSize)).^2;
            FM = std2(LVar)^2;
        else
            LVar = stdfilt(Image, ones(WSize, WSize)).^2;
            FM = stdfilt(LVar, ones(WSize, WSize)).^2;
        end
        
    case 'GLVN' % Normalized GLV (Santos97)
        if AF, FM = std2(Image)^2/mean2(Image);
        else
            U = imfilter(Image, MEANF, 'replicate');
            FM = stdfilt(Image, ones(WSize,WSize)).^2./double(U);
        end

    case 'GLVM' %Modified Graylevel variance
        U = imfilter(Image, MEANF, 'replicate');
        FM = (Image-U).^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
            
    case 'GRAE' % Energy of gradient (Subbarao92a)
        Ix = Image;
        Iy = Image;
        Iy(1:end-1,:) = diff(Image, 1, 1);
        Ix(:,1:end-1) = diff(Image, 1, 2);
        FM = Ix.^2 + Iy.^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'GRAT' % Thresholded gradient (Snatos97)
        Ix = Image;
        Iy = Image;
        Iy(1:end-1,:) = diff(Image, 1, 1);
        Ix(:,1:end-1) = diff(Image, 1, 2);
        FM = max(abs(Ix), abs(Iy));
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'GRAS' % Squared gradient (Eskicioglu95)
        Ix = Image;
        Ix(:,1:end-1) = diff(Image, 1, 2);
        FM = Ix.^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'HELM' %Helmli's mean method (Helmli2001)
        U = imfilter(Image, MEANF, 'replicate');
        R1 = U./Image;
        R1(Image==0)=1;
        index = (U>Image);
        FM = 1./R1;
        FM(index) = R1(index);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'HISE' % Histogram entropy (Krotkov86)
        if AF, FM = entropy(Image);
        else
            NHOOD = ones(WSize,WSize);
            FM = entropyfilt(Image, NHOOD);
            FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'HISR' % Histogram range (Firestone91)
        if AF, FM = max(Image(:))-min(Image(:));
        else
        NHOOD = ones(WSize,WSize);
        FM = rangefilt(Image,NHOOD);
        FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'HISV' % Variance of log Histogram (Forero2004)
        error('Not implemented yet, sorry')
        
    case 'LAPE' % Energy of laplacian (Subbarao92a)
        LAP = fspecial('laplacian');
        FM = imfilter(Image, LAP, 'replicate', 'conv');
        if AF, FM = mean2(FM.^2);
        else FM = imfilter(FM.^2, MEANF, 'replicate');
        end
                
    case 'LAPM' % Modified Laplacian (Nayar89)
        M = [-1 2 -1];
        Lx = imfilter(Image, M, 'replicate', 'conv');
        Ly = imfilter(Image, M', 'replicate', 'conv');
        FM = abs(Lx) + abs(Ly);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'LAPV' % Variance of laplacian (Pech2000)
        LAP = fspecial('laplacian');
        ILAP = imfilter(Image, LAP, 'replicate', 'conv');
        if AF, FM = std2(ILAP)^2;
        else FM = stdfilt(ILAP,ones(WSize,WSize)).^2;
        end
        
    case 'LAPD' % Diagonal laplacian (Thelen2009)
        M1 = [-1 2 -1];
        M2 = [0 0 -1;0 2 0;-1 0 0]/sqrt(2);
        M3 = [-1 0 0;0 2 0;0 0 -1]/sqrt(2);
        F1 = imfilter(Image, M1, 'replicate', 'conv');
        F2 = imfilter(Image, M2, 'replicate', 'conv');
        F3 = imfilter(Image, M3, 'replicate', 'conv');
        F4 = imfilter(Image, M1', 'replicate', 'conv');
        FM = abs(F1) + abs(F2) + abs(F3) + abs(F4);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'GRA3' %3D Laplacian (Ahmad2007)
        if AF, error('GRA3 is not an AF operator');
        end
        [M, N, L] = size(Image);
        FM = zeros(M,N,L);
        a = sqrt(2)/2; b = sqrt(3)/3;
        Ox(1,1:3,1:3) = [b a b; a 1 a; b a b];
        Ox(2,:,:) = zeros(3,3);
        Ox(3,:,:) = -Ox(1,:,:);
        Oy(1:3,1,1:3) = Ox(1,:,:);
        Oy(:,2,:) = Ox(2,:,:);
        Oy(:,3,:) = Ox(3,:,:);
        Oz(:,:,1) = Ox(1,:,:);
        Oz(:,:,2) = Ox(2,:,:);
        Oz(:,:,3) = Ox(3,:,:);
        Vx = (imfilter(Image, Ox, 'conv','symmetric'));
        Vy = (imfilter(Image, Oy, 'conv','symmetric'));
        Vz = (imfilter(Image, Oz, 'conv','symmetric'));        
        for k = 1:L
            V = sqrt(double(Vx(:,:,k)).^2+...
                double(Vy(:,:,k)).^2+double(Vz(:,:,k)).^2);
            FM(:,:,k) = imfilter(V, MEANF, 'replicate');
        end
        
    case 'LAP3' % Laplacian in 3D Window (An2008)
        if AF %error('LAP3 is not an AF operator');
            FM = mean2(focusmeasure(Image, 'LAPM', WSize));            
            return
        end
        [M,N,L] = size(Image);
        V = zeros(M, N, 3);
        FM = zeros(M, N, L);
        FM(:,:,1) = focusmeasure(im2double(Image(:,:,1)),...
            'lapm', WSize);
        for k = 2:L-1
            V(:,:,1) = focusmeasure(im2double(Image(:,:,k-1)),...
                'lapm', WSize);
            V(:,:,2) = focusmeasure(im2double(Image(:,:,k)),...
                'lapm', WSize);
            V(:,:,3) = focusmeasure(im2double(Image(:,:,k+1)),...
                'lapm', WSize);
            FM(:,:,k) = mean(V,3);
        end
        FM(:,:,L) = focusmeasure(im2double(Image(:,:,1)),...
            'lapm', WSize);
        
    case 'LBPA' %Local Binary Patterns 1 (Lorenzo2008)
        %Image = im2double(Image);
        MAPPING = getmapping(16,'riu2');
        FM = zeros(size(Image));
        FM(3:end-2,3:end-2) = lbp(Image,2,8,MAPPING,'nohist');
        %FM = lbp(Image,2,16);
        if AF, FM = 1/mean2(FM);
        else FM = 1./imfilter(FM, MEANF, 'replicate');
        end
        
    case 'SFIL' %Steerable filters (Minhas2009)
        % Angles = [0 45 90 135 180 225 270 315];
        N = floor(WSize/2);
        sig = N/3;
        [x,y] = meshgrid(-N:N, -N:N);
        G = exp(-(x.^2+y.^2)/(2*sig^2))/(2*pi*sig);
        Gx = -x.*G/(sig^2);Gx = Gx/sum(Gx(:));
        Gy = -y.*G/(sig^2);Gy = Gy/sum(Gy(:));
        R(:,:,1) = imfilter(Image, Gx, 'conv', 'replicate');
        R(:,:,2) = imfilter(Image, Gy, 'conv', 'replicate');
        R(:,:,3) = cosd(45)*R(:,:,1)+sind(45)*R(:,:,2);
        R(:,:,4) = cosd(135)*R(:,:,1)+sind(135)*R(:,:,2);
        R(:,:,5) = cosd(180)*R(:,:,1)+sind(180)*R(:,:,2);
        R(:,:,6) = cosd(225)*R(:,:,1)+sind(225)*R(:,:,2);
        R(:,:,7) = cosd(270)*R(:,:,1)+sind(270)*R(:,:,2);
        R(:,:,7) = cosd(315)*R(:,:,1)+sind(315)*R(:,:,2);
        FM = max(R,[],3);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'SFRQ' % Spatial frequency (Eskicioglu95)
        Ix = Image;
        Iy = Image;
        Ix(:,1:end-1) = diff(Image, 1, 2);
        Iy(1:end-1,:) = diff(Image, 1, 1);
        if AF, FM = mean2(sqrt(double(Iy.^2+Ix.^2)));
        else 
        RF = imfilter(Iy.^2, MEANF, 'replicate');
        CF = imfilter(Ix.^2, MEANF, 'replicate');
        FM = sqrt(double(RF) + double(CF));
        end
        
    case 'TENG'% Tenengrad (Krotkov86)
        Sx = fspecial('sobel');
        Gx = imfilter(Image, Sx, 'replicate', 'conv');
        Gy = imfilter(Image, Sx', 'replicate', 'conv');
        FM = Gx.^2 + Gy.^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'TENV' % Tenengrad variance (Pech2000)
        Sx = fspecial('sobel');
        Gx = imfilter(Image, Sx, 'replicate', 'conv');
        Gy = imfilter(Image, Sx', 'replicate', 'conv');
        G = Gx.^2 + Gy.^2;
        if AF, FM = std2(G)^2;
        else FM = stdfilt(G, ones(WSize,WSize)).^2;
        end
        
    case 'VOLA' % Vollath's correlation (Santos97)
        %Image = im2double(Image);
        I1 = Image; I1(1:end-1,:) = Image(2:end,:);
        I2 = Image; I2(1:end-2,:) = Image(3:end,:);
        if AF, 
            Image = Image.*(I1-I2);
            FM = mean2(Image);
        else
            FM = imfilter(Image.*(I1-I2), MEANF, 'replicate');            
        end
        
    case 'VOLS' % Vollath's STD (Santos97)
        %Image = im2double(Image);
        U = WSize*imfilter(Image, MEANF, 'replicate');
        I1 = Image;        
        I1(1:end-1,:) = Image(2:end,:);
        FM = Image.*I1-U.^2;
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'WAVS' %Sum of Wavelet coeffs (Yang2003)
        [C,S] = wavedec2(Image, 1, 'db6');
        H = wrcoef2('h', C, S, 'db6', 1);   
        V = wrcoef2('v', C, S, 'db6', 1);   
        D = wrcoef2('d', C, S, 'db6', 1);   
        FM = abs(H) + abs(V) + abs(D);
        if AF, FM = mean2(FM);
        else FM = imfilter(FM, MEANF, 'replicate');
        end
        
    case 'WAVV' %Variance of  Wav...(Yang2003)
        [C,S] = wavedec2(Image, 1, 'db6');
        H = abs(wrcoef2('h', C, S, 'db6', 1));
        V = abs(wrcoef2('v', C, S, 'db6', 1));
        D = abs(wrcoef2('d', C, S, 'db6', 1));
        if AF
            FM = std2(H)^2+std2(V)+std2(D);
        else
            NHOOD = ones(WSize, WSize);
            FM = stdfilt(H,NHOOD).^2+stdfilt(V,NHOOD).^2+...
                stdfilt(D, NHOOD).^2;
        end
        
    case 'WAVR' %Ratio of wavelet coefficients (Xie2006)
        [C,S] = wavedec2(Image, 3, 'db6');
        H = abs(wrcoef2('h', C, S, 'db6', 1));   
        V = abs(wrcoef2('v', C, S, 'db6', 1));   
        D = abs(wrcoef2('d', C, S, 'db6', 1)); 
        A1 = abs(wrcoef2('a', C, S, 'db6', 1));
        A2 = abs(wrcoef2('a', C, S, 'db6', 2));
        A3 = abs(wrcoef2('a', C, S, 'db6', 3));
        A = A1 + A2 + A3;
        WH = H.^2 + V.^2 + D.^2;
        if AF, 
            WH = mean2(WH);
            WL = mean2(A);
            FM = WH/WL;
        else
            
            WH = imfilter(WH, MEANF, 'replicate');
            WL = imfilter(A.^2, MEANF, 'replicate');
            FM = WH./WL;
        end
    case 'WAVC' %Curvelet transform [Minhas2011]
        C = fdct_wrapping(Image, 1, 1, 2, 8);        
        C1 = imresize(C{1}{1},size(Image));
        C2 = imresize(C{2}{1},size(Image))+...
            imresize(C{2}{2},size(Image))+...
            imresize(C{2}{3},size(Image))+...
            imresize(C{2}{4},size(Image))+...
            imresize(C{2}{4},size(Image))+...
            imresize(C{2}{5},size(Image))+...
            imresize(C{2}{6},size(Image))+...
            imresize(C{2}{7},size(Image))+...
            imresize(C{2}{8},size(Image));        
        if AF, FM = sum(C2(:)./C1(:));
        else FM = imfilter(C2./C1,MEANF,'replicate');
        end
    otherwise
        error('Unknown measure %s',upper(Measure))
end
 end
%************************************************************************
function fm = AcMomentum(Image)
[M, N] = size(Image);
Hist = imhist(Image)/(M*N);
Hist = abs((0:255)-255*mean2(Image))'.*Hist;
fm = sum(Hist);
end
%******************************************************************
function F = TchebiFocus(Image)
Ord = 2;                    %Order of Tchebichef polynomials
[M, N] = size(Image);
Image = Image./sqrt(M*N*mean2(Image.^2));
I = Image - mean2(Image);
P1 = TchebichefPoly(M, Ord); % P1 = Tchebichef polynomials of x-axis
P2 = TchebichefPoly(N, Ord); % P2 = Tchebichef polynomials of y-axis
M2 = size(P1, 1);            % Determine size of P1
N2 = size(P1, 1);            % Determine size of P2
B = P1*I*P2';                % Compute Tchebichef moments
B = fliplr(triu(fliplr(B),0)); % extract low order Tchebichef moments
I2 = I.^2;
B2 = B(1:M2,1:N2).^2;
F = abs( sum(I2(:)) - sum(B2(:)))/abs(sum(B2(:)));
end
%*******************************************************************
function F = TchebichefPoly(N, Ord)

Pnx = double(zeros(Ord+1,N));           % Tchebichef polynomials
for x = 0:N-1
    x1 = x + 1;
    Pnx(0+1,x1) = 1.0;                  % P0 = 1.0
    Pnx(1+1,x1) = (1 - N + 2.0*x) / N;  % P1 = (1-N+2x); divided by N to map Tchebichef polynomials within [-1,1]
    for p = 2:Ord,
        p1 = p + 1;
        Pnx(p1,x1) = ( double(2*p-1)*Pnx(1+1,x1)*Pnx(p1-1,x1) - ...
            double((p-1))*(1.0-double((p-1)*(p-1))/double(N*N))*Pnx(p1-2,x1) ) / double(p);
    end
end

R = Rho(N);
for x = 0:N-1
    x1 = x + 1;
    for p = 0:Ord
        p1 = p + 1;
        Pnx(p1,x1) = Pnx(p1,x1) / sqrt(R(p1,1));
    end
end
F = Pnx;
end
%*********************************************************************
function R = Rho(N)

% This function is created to compute the rho (normalized factor) 
% of Tchebichef polynomials.

R = zeros(N,1);
R(0+1) = 1.0;

for p = 1:N-1
    p1 = p + 1;
    R(p1,1) = R(p1-1,1) * (1.0 - double((p*p))/double((N*N)));  % (1-p^2/N^2)
end
for p = 0:N-1
    p1 = p + 1;
    R(p1,1) = R(p1,1) * double(N) / double((2*p + 1));      % N(1-P^2/N^2)/(2p+1)
end
end
%***********************************************************************
function fm = DctRatio(M)
MT = dct2(M).^2;
fm = (sum(MT(:))-MT(1,1))/MT(1,1);
end
%************************************************************************
function fm = ReRatio(M)
M = dct2(M);
fm = (M(1,2)^2+M(1,3)^2+M(2,1)^2+M(2,2)^2+M(3,1)^2)/(M(1,1)^2);
end

%******************************************************************
function FM = EigenFocus(Image)

K = 3; %Method parameter
[M, N] = size(Image);
Image = Image./sqrt(M*N*mean2(Image.^2));
I = Image - mean2(Image);
covariance = (I*I')/(N*M-1);
[~, V] = eig(covariance);
var = diag(V);
[~, rindices] = sort(-var);
var = var(rindices);
FM = sum(var(1:K));
end
%******************************************************************
