k=1;
for i=80:130
    for j=40:60
        v(k)= output_2d_s(i,j);
        k = k+1;
    end
end
sb = mean(v(:))
sigb = std(v(:))
snr=sb/sigb