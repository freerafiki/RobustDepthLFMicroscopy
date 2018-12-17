function [fs, aif] = read_fs_aif(path)

fprintf('\nReading the focal stack and the central image..\n');
% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
cd '/Users/Palma/Documents/Valencia/Images/Type1/FocalStack3/';
cd 'fs';
imagefiles = dir('*.bmp');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   [strings, f_part] = strsplit(currentfilename, 'f=');
   [strings2, ext] = strsplit(strings{2}, '.bmp');
   fl_string = strings2{1};
   fl = str2num(fl_string);
   fs{ii,1} = rgb2gray(currentimage);
   fs{ii,2} = fl;
end

cd '/Users/Palma/Documents/Valencia/Images/Type1/FocalStack3/';
cd 'aif';
imagefiles = dir('*.bmp');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   aif{ii} = currentimage;
end
if size(aif,2) == 1
    aif = aif{1};
end

cd '/Users/Palma/Documents/Valencia/Code/matlab_depth/';
fprintf('Done!\n');