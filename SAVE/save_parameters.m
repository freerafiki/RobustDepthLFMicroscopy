function save_parameters(folder_path, img_type, pitch, ini, fin, a, C0, offset, step_pix, ...
    path, superpixels_size, window_size, alpha_DFD, alpha_DFC, matting, gamma_DFD_s1, ...
    gamma_DFD_s2, beta_DFD, gamma_DFC_s1, gamma_DFC_s2, beta_DFC, phi, sp_type, wsPP)

nameFile = 'parameters.txt';
fileID = fopen(strcat(folder_path, nameFile),'w');
fprintf(fileID,'PARAMETERS\n\n');
fprintf(fileID,'Imgtype=%d', img_type);
% TYPE 1 - Fibers
% TYPE 3 - Bicho
% TYPE 4 - Pieces
% TYPE 5 - Telescope
% TYPE 6 - Test
% TYPE 7 - Flies
% TYPE 8 - Chip
% TYPE 9 - Flower (with matting!) - Experimental
if img_type == 1
    fprintf(fileID,' fibers');
elseif img_type == 2
    fprintf(fileID,' unknown');
elseif img_type == 3
    fprintf(fileID,' bicho');
elseif img_type == 4
    fprintf(fileID,' mcf');
elseif img_type == 5
    fprintf(fileID,' telescope');
elseif img_type == 6
    fprintf(fileID,' test');
elseif img_type == 7
    fprintf(fileID,' flies');
elseif img_type == 8
    fprintf(fileID,' chips');
elseif img_type == 9
    fprintf(fileID,' flower (with matting) - experimental');
end
fprintf(fileID,'\npitch=%d, ini=%d, fin=%d, a=[%d %d %d], C0=[%d, %d]\n', pitch, ini, fin, a(1), a(2), a(3), C0(1), C0(2));
fprintf(fileID,'offset=%5.1f, step_pix=%3.1f\n', offset, step_pix);
fprintf(fileID,'Superpixel technique=');
if sp_type == 1
    fprintf(fileID,'vl_slic from MATLAB, superpixels size = %d\n',superpixels_size);
else
    fprintf(fileID,'slicmex from EPFL, number of sp = 500, compactness = 20\n');
end
fprintf(fileID,'window_size=%d, merging_parameter=%1.3f, matting=%d, window_size_post-proc=%d\n', window_size, phi, matting, wsPP);
fprintf(fileID,'\nDEFOCUS\n');
fprintf(fileID,'alphaDFD=%1.3f, gammaDFD_s1=%1.3f, gammaDFD_s2=%1.3f, betaDFD=%1.3f\n', alpha_DFD, gamma_DFD_s1, gamma_DFD_s2, beta_DFD);
fprintf(fileID,'\nCORRESPONDENCES\n');
fprintf(fileID,'alphaDFC=%1.3f, gammaDFC_s1=%1.3f, gammaDFC_s2=%1.3f, betaDFC=%1.3f\n', alpha_DFC, gamma_DFC_s1, gamma_DFC_s2, beta_DFC);
fprintf(fileID,'\nCalculated on image at ');
fprintf(fileID,path);
