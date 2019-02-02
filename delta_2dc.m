function del_2d=delta_2d(i,j,da,dl,g,g_prime)
% i is the row number and j is the A-lines and da and dl are axial and
% lateral displacements
del_2d=abs(g(i,j)-g_prime(i+da,j+dl));
end