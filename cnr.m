
k=1;
for i=180:230
    for j=15:25
        v(k)= output_2d_s(i,j);
        k = k+1;
    end
end
sb = mean(v(:))
sigb = var(v(:))
q=1;
for i=160:220
    for j=45:55
        z(q)= output_2d_s(i,j);
        q = q+1;
    end
end
st = mean(z(:))
sigt = var(z(:))

Cnr = sqrt((2*(sb - st)^2)/(sigb^2 + sigt^2));
Cnr

                