function s_2d=S_2d(da_i,dl_i,da_prev,dl_prev)
% da_i is the axial displacement at the ith row-pixel and dl_i is the lateral one
% usually, they're a scalar.
% da_prev and dl_prev are the axial and lateral displacement of the previous row
% let's say we want to calculate 10th row so we need value at the 9th
% row;both axial and lateral.
s_2d=(da_i-da_prev).^2+(dl_i-dl_prev).^2;
end