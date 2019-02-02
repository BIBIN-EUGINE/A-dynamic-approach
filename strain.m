function [strain_img]=strain(img,m)
%img is displacement matrix and m is the size of filter.
%resize to calculate faster
img=imresize(img,0.2);

if(mod(m,2)==0)
    return
end
n=(m-1)/2; %how many rows we should skip at the end and the begining.
[M,N]=size(img);

%img_vec=reshape(img, M*N,1);
strain_img=zeros(size(img));
strain_col=zeros(M,1);
idx=[1:m]';
i=M*N-n+1
for j=1:N
    column=img(:,j);
    for i=n+1:M-n
        tmp=polyfit(idx,column(i-n:i+n),1);
        strain_col(i)=tmp(1);
    end
    strain_img(:,j)=strain_col;
end

%if you want to only use absolute values of the slope, uncomment next line.
strain_img=abs(strain_img);

%normalization
m=mean(strain_img(:));
s=std(strain_img(:));

max_limit=m+3*s;
min_limit=m-3*s;

strain_img(strain_img>max_limit)=max_limit;
strain_img(strain_img<min_limit)=min_limit;

%reshape
strain_img=reshape(strain_img, M, N);
%SNR=m/s;
return
end