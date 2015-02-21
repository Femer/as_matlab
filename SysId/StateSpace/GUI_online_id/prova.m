fileID = fopen('test.txt');

C = textscan(fileID, '%s %s %s %s', 1);
celldisp(C)

C = textscan(fileID, '%f %f %f %f');
celldisp(C)

fclose(fileID);
