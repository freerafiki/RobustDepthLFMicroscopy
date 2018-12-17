function [z,R] = sff_v2(fs, varargin)
% Shape from focus (SFF)
% SINTAX:
%
%   Z  = sff(IMAGES)   
%   Z  = sff(IMAGES, option1, value1, option2, ...)
%   [Z,R] = sff(IMAGES, ...);
%
% DESCRIPTION:
% Estimate depthmap and reliability measure from defocused image 
% sequence using shape-from-focus (SFF). According to [2], depth for 
% pixels with R<20 dB is not reliable and should be descarded.
% 
% INPUTS:
% IMLIST,     is a 1xP cell a array where each cell is a 
%             string with the image path  corresponding to 
%             one frame of the focus sequence.
%
% Options and their values may be any of the following 
%(default value in parenthesis):
% 'fmeasure',   A string with the focus measure operator. (GLVM).
% 'filter',     A scalar with the size of median filter. (0).
% 'interp',     A logic flag to turn Gaussian interpolation on/off (true)
% 'fsize',      An integer with the size of the focus mesure window. (9)
% 'focus',      A vector with the focus position of each image. (1:P)
% 'ROI',        Scene ROI as a rectangle [xo yo W H]. Default is [].
%
% OUTPUTS:
% Z,        is a MxN matrix with the reconstructed depthmap obtained using
%           SFF as originally proposed in [1]
% R,        is a MxN matrix with the reliability measure (in dB) of the
%           estimated depth computed as proposed in [2]
%
% References:
% [1] S.K. Nayar and Y. Nakagawa, PAMI, 16(8):824-831, 1994. 
% [2] S. Pertuz, D. Puig, M. A. Garcia, Reliability measure for shape-from
% focus, IMAVIS, 31:725�734, 2013.
%
%S. Pertuz
%Jan/2016.

for j = 1:size(fs,4)
    fs(:,:,:,j) = uint8(fs(:,:,:,j));
end

opts.RGB = 0;
opts.INTERP = 1;
opts.fmeasure = 'LAPD';
opts.filter = 0;
opts.nhsize = 9;
images = zeros(size(fs,1), size(fs, 2), size(fs, 4));
for j = 1:size(fs,4)
    images(:,:,j) = rgb2gray(uint8(fs(:,:,:,j)));
    opts.focus(j) = j;
end
opts.ROI = [];
opts.size = [size(fs,1), size(fs,2), size(fs,4)];


% Stack size:
M = opts.size(1);
N = opts.size(2);
P = opts.size(3);
tic
% Read images:
%{
tic
fprintf('Reading      ')
images = zeros(M, N, P,'uint8');
for p = 1:P
    im = imread(imlist{p});
    
    if ~isempty(opts.ROI), im = imcrop(im, opts.ROI);
    end
    if opts.RGB, im = rgb2gray(im);
    end    
    images(:,:,p) = im;
    fprintf('\b\b\b\b\b[%2.2i%%]',round(100*p/P))
end
clear imdata
t1 = toc;
%}
t1 = toc;
%Compute focus measure volume:
tic
fprintf('\nFmeasure      ')
fm = zeros(M, N, P); 

for p = 1:P
    fm(:,:,p) = focusmeasure(im2double(images(:,:,p)),...
        opts.fmeasure, opts.nhsize);
    fprintf('\b\b\b\b\b[%2.2i%%]',round(100*p/P))
end
t2 = toc;

% Estimate depthmap
fprintf('\nDepthmap ')
if opts.INTERP&&(nargout==2) 
    [I, zi, s, A] = gauss3P(opts.focus, fm);
    z = zi;
    z(z>max(opts.focus)) = max(opts.focus);
    z(z<min(opts.focus)) = min(opts.focus);
elseif (nargout==2)
    [I, zi, s, A] = gauss3P(opts.focus, fm);    
    z = opts.focus(I);
else
    [~, I] = max(fm,[ ], 3);
    z = opts.focus(I);
end

fmax = opts.focus(I);
z(isnan(z)) = fmax(isnan(z));
t3 = toc;
fprintf('[100%%]\n' )


% Median filter:
if opts.filter~=0
    fprintf('Smoothing ')
    z = medfilt2(z, [opts.filter opts.Filter]);
    fprintf('[100%%]\n')
end

tic
% Reliability measure
if nargout==2
    fprintf('Rmeasure      ')
    err = zeros(M, N);
    
    %Compute fitting error:
    for p = 1:P
        err = err + abs( fm(:,:,p) - ...
            A.*exp(-(opts.focus(p)- zi).^2./(2*s.^2)));
        fprintf('\b\b\b\b\b[%02d%%]',round(100*p/P))
    end
    
    h = fspecial('average', opts.nhsize);
    err = imfilter(err, h);
    
    R = 20*log10(P*fmax./err);
    mask = isnan(zi);
    R(mask|R<0|isnan(R)) = 0;
    fprintf('\n')
end
t4 = toc;
% Print output
tt = t1 + t2 + t3 + t4;
fprintf('\nElapsep time: %1.2f s:\n',tt)
fprintf('Reading (%1.1f%%).\n',100*t1/tt);
fprintf('Fmeasure (%1.1f%%).\n',100*t2/tt);
fprintf('Depthmap (%1.1f%%).\n',100*t3/tt);

if nargout==2
    fprintf('Rmeasure (%1.1f%%).\n',100*t4/tt);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I, u, s, A] = gauss3P(x, Y)
% Closed-form solution for Gaussian interpolation using 3 points
% Internal parameter:
STEP = 2;
%%%%%%%%%%%%%%%%%%%
[M,N,P] = size(Y);
[~, I] = max(Y,[ ], 3);
[IN,IM] = meshgrid(1:N,1:M);
Ic = I(:);
Ic(Ic<=STEP)=STEP+1;
Ic(Ic>=P-STEP)=P-STEP;
Index1 = sub2ind([M,N,P], IM(:), IN(:), Ic-STEP);
Index2 = sub2ind([M,N,P], IM(:), IN(:), Ic);
Index3 = sub2ind([M,N,P], IM(:), IN(:), Ic+STEP);
% Index1(I(:)<=STEP) = Index3(I(:)<=STEP);
% Index3(I(:)>=STEP) = Index1(I(:)>=STEP);
x1 = reshape(x(Ic(:)-STEP),M,N);
x2 = reshape(x(Ic(:)),M,N);
x3 = reshape(x(Ic(:)+STEP),M,N);
y1 = reshape(log(Y(Index1)),M,N);
y2 = reshape(log(Y(Index2)),M,N);
y3 = reshape(log(Y(Index3)),M,N);

c = ( (y1-y2).*(x2-x3)-(y2-y3).*(x1-x2) )./...
    ( (x1.^2-x2.^2).*(x2-x3)-(x2.^2-x3.^2).*(x1-x2) );
b = ( (y2-y3)-c.*(x2-x3).*(x2+x3) )./(x2-x3);
a = y1 - b.*x1 - c.*x1.^2;

s = sqrt(-1./(2*c));
u = b.*s.^2;
A = exp(a + u.^2./(2*s.^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Opts = parseInputs(imlist, varargin)                                   %

P = length(imlist);

input_data = inputParser;
input_data.CaseSensitive = false;
input_data.StructExpand = true;

input_data.addOptional('filter', 0, @(x) isnumeric(x)&& isscalar(x));

input_data.addOptional('fmeasure', 'glvm', @(x) ischar(x));

input_data.addOptional('focus', 1:P, @(x) isnumeric(x));

input_data.addOptional('interp', true, @(x) islogical(x) &&isscalar(x));

input_data.addOptional('nhsize', 9, @(x) isnumeric(x) && isscalar(x));

input_data.addOptional('ROI', [], @(x) isnumeric(x) && numel(x)==4);

parse(input_data, varargin{:});

im = imread(imlist{1});
Opts.RGB = (numel(size(im))==3);
Opts.INTERP = input_data.Results.interp;
Opts.fmeasure = input_data.Results.fmeasure;
Opts.filter = round(input_data.Results.filter);
Opts.nhsize = round(input_data.Results.nhsize);
Opts.focus = input_data.Results.focus;
Opts.ROI = input_data.Results.ROI;
Opts.size = [size(im), length(imlist)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    