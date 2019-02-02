clear all
% close all

%pre compression image
load 'rf01.mat'
Im1 = RfDataDouble(1:1700,:);
maxIm = max(Im1(:));
Im1 = Im1/maxIm;
figure(1);
imagesc(Im1);
colormap gray

%post compression image 
load 'rf03.mat'
Im2 = RfDataDouble(1:1700,:);
Im2 = Im2/maxIm;
figure(2);
imagesc(Im2);
colormap gray

%clear signals_matrix
% Im1 = [1 2 3;4 5 6;7 8 9;10 3 5;3 7 8];
% Im2 = [1 2 3;7 8 9;4 5 6;3 7 8;10 3 5];
m = length(Im1(:,1));%taking the first column length
n = length(Im1(1,:));%Taking the first row length

dmax = 0;
neg_dmin = 100;
alpha = 0.15; % weight factor
a = [0 0 0]; % matrix for finding the minima
for i = 1:n
    g1(:,1) = Im1(:,i);%Pre compression image pixel value storing
    g2(:,1) = Im2(:,i);%Post compression image pixel value storing
    disp(i)
%     g1(:,1) = [1 2 3;4 5 6;7 8 9];
%     g2(:,1) = [1 2 3;7 8 9;4 5 6];
%     m = length(g1);
%     n = length(g2);
    c = NaN(m,neg_dmin); % cost function
    c(1,neg_dmin) = abs(g1(1) - g2(1));
    for j = 2:neg_dmin-1
        for k = neg_dmin-(j-1):neg_dmin
            
            a(1,1) = c(j-1,k-1) + alpha * (((-(neg_dmin-k))-(-(neg_dmin-k-1)))^2);
            a(1,2) = c(j-1,k) + alpha * (((-(neg_dmin-k)) - (-(neg_dmin-k)))^2);
            if (k== neg_dmin)
                a(1,3) = 1000000;
            else
                a(1,3) = c(j-1,k+1) + alpha * (((-(neg_dmin-k)) - (-(neg_dmin-k+1)))^2);
            end
            [minimum,index_1] = min(a);
            if(index_1 == 1)
                M(j,k) = k-1;
            elseif(index_1 == 2)
                M(j,k) = k;
            else
                M(j,k) = k+1;
            end
            delta = deltafn(j,neg_dmin,k,g1(:,1),g2(:,1));
            c(j,k) = minimum + delta; % cost function
        end
    end
            
    for j = neg_dmin:m-1
        for k = 1:neg_dmin
            if(k==1)
                a(1,1) = 100000;
            else
                a(1,1) = c(j-1,k-1) + alpha * (((-(neg_dmin-k))-(-(neg_dmin-k-1)))^2);
            end
            a(1,2) = c(j-1,k) + alpha * (((-(neg_dmin-k)) - (-(neg_dmin-k)))^2);
            if (k== neg_dmin)
                a(1,3) = 1000000;
            else
                a(1,3) = c(j-1,k+1) + alpha * (((-(neg_dmin-k)) - (-(neg_dmin-k+1)))^2);
            end
            [minimum,index_1] = min(a);
            if(index_1 == 1)
                M(j,k) = k-1;
            elseif(index_1 == 2)
                M(j,k) = k;
            else
                M(j,k) = k+1;
            end
            c(j,k) = minimum + deltafn(j,neg_dmin,k,g1(:,1),g2(:,1)); % cost function
        end
         if(j == m-1)
            [Q, last_index]   = min(c(j,:));
         end
    end 
    d(m-1,i) = last_index;
    for j = m-2:-1:1
        d(j,i) = M(j+1,d(j+1,i)); % displacement map
    end  
end
figure(4),imagesc(d),title('Axial Displacement'),colormap gray,colorbar();
strain_grad = gradient(d);
figure,imagesc(strain_grad),colormap gray,colorbar();
filtered_displacement = medfilt2(d,[3 3]);
figure(5),imagesc(filtered_displacement),title('Axial Displacement after filtering'),colormap gray,colorbar();
strain_image=strain(filtered_displacement,5);
figure(6), imagesc(strain_image),title('Strain Image'),colormap gray,colorbar();








 