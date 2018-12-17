function save_final_depth(folder_path, depth, name)

name1 = strcat(folder_path, name, '.png');  
name2 = strcat(folder_path, name, '_no_norm.png'); 
maxval = max(max(depth));
imgtosave = depth./maxval;
imwrite(imgtosave,name1);
imwrite(uint8(depth),name2);