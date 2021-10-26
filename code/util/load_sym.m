function sym = load_sym(symdata)

fileID = fopen(symdata, 'r');
formatSpec = '%f';
sym = fscanf(fileID, formatSpec);

end