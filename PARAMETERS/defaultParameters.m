function [pitch, C0, a, ini, fin, offset, step_pix] = defaultParameters(type)

% TYPE1 = microscope - fibers
if type == 1 
    pitch=683;%telescope444; %473; % 682
    pt=floor(pitch/2);
    C0=[478, 831]; %2 %%   [472, 848]; %3 %%    %;%[295, 536];%telescope[344, 496];
    a=[2 3 2];
    ini = 666;
    fin = 697;
    offset=70; %%%focus offset of the optical system
    step_pix=13.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel
% BICHO
elseif type == 3
    pitch=473;
    pt=floor(pitch/2);
    C0=[288, 514]; % [287, 513];%[295, 536];%telescope[344, 496];
    a=[2 3 2];
    ini = 445;
    fin = 465;
    offset=70; %%%focus offset of the optical system
    step_pix=13.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel   
% stranger things
elseif type == 4
    pitch=687;
    pt=floor(pitch/2);
    C0=[486, 851]; % [287, 513];%[295, 536];%telescope[344, 496];
    a=[2 3 2];
    ini = 680;
    fin = 699;
    offset=70; %%%focus offset of the optical system
    step_pix=13.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel  
% TYPE5 = telescope
elseif type == 5 
    pitch=483; % 447;%telescope444; %473; % 682
    pt=floor(pitch/2);
    C0=[783, 1104]; % [350, 500];%[295, 536];%telescope[344, 496];
    a=[2 3 2];
    ini = 493;
    fin = 507;
    offset=200; %%%focus offset of the optical system
    step_pix=3.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel
% TYPE6 = test
elseif type == 6
    pitch=473;
    pt = floor(pitch/2);
    C0 = [343, 567];
    a=[2 3 2];
    ini = 448;
    fin = 456;
    offset = 70;
    step_pix = 13.5;
% Flies
elseif type == 7
    pitch=487;
    pt = floor(pitch/2);
    C0 = [397, 580];
    a=[2 3 2];
    ini = 505;
    fin = 515;
    offset = 70;
    step_pix = 13.5;
% Chip
elseif type == 8
    pitch=445;
    pt = floor(pitch/2);
    C0 = [307 997]; %%%[294 596]; % chip21 ------ %%% [449, 585];% chip 16 %%%
    a=[2 3 2];
    ini = 445;
    fin = 473;
    offset = 0; %70;
    step_pix = 1; % 13.5;
% flower with matting
elseif type == 9
    pitch=455;
    pt = floor(pitch/2);
    C0 = [397, 580];
    a=[2 3 2];
    ini = 445;
    fin = 460;
    offset = 0; %70;
    step_pix = 1; % 13.5;
% SYNTHETIC
elseif type == 10
    pitch=683;%telescope444; %473; % 682
    pt=floor(pitch/2);
    C0=[478, 831]; %2 %%   [472, 848]; %3 %%    %;%[295, 536];%telescope[344, 496];
    a=[2 3 2];
    ini = 683;
    fin = 711;
    offset=70; %%%focus offset of the optical system
    step_pix=13.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel
% ZEBRAFISH
elseif type == 11
    pitch=445;
    pt = floor(pitch/2);
    C0 = [307 997]; 
    a=[2 3 2];
    ini = 445;
    fin = 468;
    offset = 0; %70;
    step_pix = 1; % 13.5;
else
    pitch=447; % 682
    pt=(pitch/2);
    %center of EI(1)
    C0=[295, 536];  %478, 831	
    %EI in each row
    a=[2 3 2];
    %refocusing ini and fin
    ini=436; %670;
    fin=456; %694;
    offset=70; %%%focus offset of the optical system
    step_pix=13.5; %14.5; % step_pix=27.5/scale; % 5-4-5-4 structure is 27.5 micrometer/pixel
end
