% *************************************************************************
% v3 some speed issues and strain image calculation integration
% v4 includes lateral displacement instead of only axial.
% v5 resolves speed issues
% *************************************************************************

clear all
% close all

%inputs
w=.15; %this is the regularization weight
% g=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15];
% g_prime=[0;2;3;4;6;7;8;9;10;9;10;11;12;13;14];
load 'rf01.mat'
Im1 = RfDataDouble(1:1700,:);
Im1=Im1/max(Im1(:));
% Im1=[1;2;3;4;5;6;7;8;9;10];
%Im1=imresize(Im1,0.5);

load 'rf03.mat'
Im2 = RfDataDouble(1:1700,:);
Im2=Im2/max(Im2(:));
% Im2=[1;1;2;3;4;5;6;7;8;9];
%Im2=imresize(Im2,0.5);


%has to include zero!!!!!!!
da_max=1;
da_min=-100;
da=[da_min:1:da_max]; %axial row vector
dl_max=4;
dl_min=-1;
dl=[dl_min:1:dl_max];%lateral row vector

% axial displacement variables
len_disp_a=size(da,2);
neg_disp_a=size(find(da<0),2);
positive_disp_a=size(find(da>0),2);
zero_da_idx=find(da==0);

% lateral displacement variables
len_disp_l=size(dl,2);
neg_disp_l=size(find(dl<0),2);
positive_disp_l=size(find(dl>0),2);
zero_dl_idx=find(dl==0);


% ensure image is bigger than dmax and dmin.
last_row=size(Im1,1);
num_cols=size(Im1,2);
num_disp_a=size(da,2); 
num_disp_l=size(dl,2);
disp_array_a=[1:1:num_disp_a];
disp_array_l=[1:1:num_disp_l];

output_2d_a=zeros(last_row,num_cols);
output_2d_l=zeros(last_row,num_cols);
disp 'calculating displacement map'

%not efficient, but easiest to test:
g=Im1;
g_prime=Im2;

%initializations
C = NaN(last_row, num_cols, num_disp_a, num_disp_l);
a = NaN(last_row, num_disp_a, num_disp_l);
M = NaN(last_row, num_disp_a, num_disp_l,2); %1 is for del_a, 2 is for del_l
disp_candidates=[repmat([1:num_disp_a]',num_disp_l,1) reshape(repmat([1:num_disp_l]',1,num_disp_a)',num_disp_l*num_disp_a,1)];


%for all columns:
first_valid_dl_idx=zero_dl_idx+1; %+1 because first column is also calculated.
last_valid_dl_idx=num_disp_l;

%iterate over all kth columns on the image
for k=1:num_cols
    k=k
    
    %control lateral displacement
    if(k<=neg_disp_l+1)
        old_first_valid_dl_idx=first_valid_dl_idx;
        first_valid_dl_idx=first_valid_dl_idx-1;
    end
    if(k>num_cols-positive_disp_l)
        old_last_valid_dl_idx=last_valid_dl_idx;
        last_valid_dl_idx=last_valid_dl_idx-1;
    end
    
    %Caculate first row. It doesn't use the previous di-1
    tmp=delta_2dc(1,k,da(da>=0), dl(first_valid_dl_idx:last_valid_dl_idx), g, g_prime);
    calc_range_a=[neg_disp_a+1:num_disp_a];
    calc_range_l=[first_valid_dl_idx:last_valid_dl_idx];
    C(1,k,calc_range_a, calc_range_l)=tmp;
    
    %for all other rows:
    first_valid_da_idx=zero_da_idx;
    last_valid_da_idx=num_disp_a;
    
    
    %for each row.
    for i=2:last_row
        %lateral displacement restrictions are already accounted for,
        %outside the loop.
        %if the row have restrictions in the axial displacement evaluation
        if(i<=neg_disp_a+1)
            old_first_valid_da_idx=first_valid_da_idx;
            first_valid_da_idx=first_valid_da_idx-1;
        end
        if(i>last_row-positive_disp_a)
            old_last_valid_da_idx=last_valid_da_idx;
            last_valid_da_idx=last_valid_da_idx-1;
        end
        
        %calculate cost for each axial and lateral displacement of this row/column
        for r=first_valid_dl_idx:last_valid_dl_idx %lateral displ.
            for j=first_valid_da_idx:last_valid_da_idx %axial displ.
                tmp=NaN(num_disp_a,num_disp_l);
                tmp2=S_2dc(da(j), dl(r), da(disp_candidates(:,1)), dl(disp_candidates(:,2))); %all values of da
                tmp=reshape(C(i-1,k,:,:), num_disp_a*num_disp_l, 1) + w*tmp2';
                                
                [a(i,j,r),idx]=min(tmp(:));
                
                del_a=disp_candidates(idx,1);
                del_l=disp_candidates(idx,2);
                M(i,j,r,:)= [del_a,del_l];
                C(i,k,j,r)=delta_2dc(i,k,da(j),dl(r),g,g_prime)+a(i,j,r);
            end
        end
        
    end

    %Traceback with indexes
    disp_m_of_row_col=reshape(C(last_row,k,:,:),num_disp_a,num_disp_l);
    [row_idx,col_dx]=find(disp_m_of_row_col==min(disp_m_of_row_col(:)),1,'first');
    D_idx=zeros(last_row,2); %1st column in axial and second is lateral
    D_idx(last_row,:)=[row_idx,col_dx];
    for i=1:last_row-1
        j=last_row-i;
        D_idx(j,:)=M(j+1,D_idx(j+1,1), D_idx(j+1,2), :);
    end

    Da=da(D_idx(:,1))';
    Dl=dl(D_idx(:,2))';
    output_2d_a(:,k)=Da;    
    output_2d_l(:,k)=Dl;
end

figure;imagesc(output_2d_a); colormap gray; colorbar; title('Axial Displacement');
output_2d_a_median=medfilt2(output_2d_a);
figure;imagesc(output_2d_a_median); colormap gray; colorbar; title('Axial Displacement After Filtering');

figure;imagesc(output_2d_l); colormap gray; colorbar; title('Lateral Displacement');

disp 'calculating strain. wait a min.'
output_2d_s=strain(output_2d_a,5);
figure; imagesc(output_2d_s); colormap hot; colorbar;title('Axial Strain Image');
output_2d_s=strain(output_2d_l,5);
figure; imagesc(output_2d_s); colormap hot; colorbar;title('Lateral Strain Image');




