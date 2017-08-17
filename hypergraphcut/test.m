f = fopen('Flori-mat.txt');
formatSpec = '%f';
A = fscanf(f, formatSpec);
disp(length(A));
INC = reshape(A, 491190, 128);
INC = sparse(INC)
save('Flori.mat', 'INC');
w = ones(491190,1);
save('Flori.mat', 'w', '-append');