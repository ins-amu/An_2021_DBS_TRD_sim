% connectivity matrix generation (weight matrix, tract-length matrix calculation)
% based on 100k streamlines (node-to-node connectivity) 
% after running the code below, save the variables wei_st_save_11, leng_st_save_11 as .mat files in each patient folder

clear all

% load fiber-tracts streamlines from TVB reconstruction output
% (TVB recon pipeline reconstructs streamlines (coordinates) based on DTI space. Thus, these should be converted into T1 space using the Flirt command in FSL, then loaded.)

for i1=0:1:99999
    file_n=sprintf('patient folder path.../SL_100k_T1_space/stream_lines-%d.txt',i1);
    data{i1+1}=dlmread(file_n,'',1,0);
end    
stream_lines=data;

% load cortical/subcortical surfaces from surf_const.m output

tri_c2=load('patient folder path.../tri_ig0_cort.mat'); tri_c3=tri_c2.tri_c2;
tri_s2=load('patient folder path.../tri_ig0_subcort.mat'); tri_s3=tri_s2.tri_s2;
ver_c2=load('patient folder path.../ver_ig0_cort.mat'); ver_c3=ver_c2.ver_c2;
ver_s2=load('patient folder path.../ver_ig0_subcort.mat'); ver_s3=ver_s2.ver_s2;
[nf_c,nv_c] = reducepatch(tri_c3,ver_c3,1/30);
[nf_s,nv_s] = reducepatch(tri_s3,ver_s3,1/30);
ver_all=[nv_c;nv_s];
tri_all=[nf_c;nf_s+length(nv_c)];

% weight/tract-length matrix generation

wei_st=sparse(zeros(size(ver_all,1),size(ver_all,1)));
leng_st=sparse(zeros(size(ver_all,1),size(ver_all,1)));
wei_final=sparse(zeros(size(ver_all,1),size(ver_all,1)));
leng_final=sparse(zeros(size(ver_all,1),size(ver_all,1)));

for num1=1:length(stream_lines)
        
    wei_st1=sparse(zeros(size(ver_all,1),size(ver_all,1)));
    wei_st2=sparse(zeros(size(ver_all,1),size(ver_all,1)));
    str_coor=stream_lines{num1};
    st_coor=str_coor(1,:); end_coor=str_coor(end,:);
    st_dis=sqrt(sum((ver_all-st_coor).^2,2)); end_dis=sqrt(sum((ver_all-end_coor).^2,2));
    st_ind=find(st_dis<10); end_ind=find(end_dis<10); % find vertices located within 10mm of both ends of each streamline
    
    if isempty(st_ind) || isempty(end_ind)
%         [num1,0]
    else
%         [num1,length(st_ind),length(end_ind)]
        wei_st1(st_ind,end_ind)=1;
        wei_st2=wei_st1+wei_st1'; wei_st2(wei_st2>1)=1;
        leng_st1=wei_st2*size(str_coor,1);

        wei_st=wei_st+wei_st2;
        leng_st=leng_st+leng_st1;
    end    
    
    if mod(num1,10000)==0 % to speed up
        wei_final=wei_final+wei_st;
        leng_final=leng_final+leng_st;
        wei_st=sparse(zeros(size(ver_all,1),size(ver_all,1)));
        leng_st=sparse(zeros(size(ver_all,1),size(ver_all,1)));
    end
    
end

% tract-length matrix
leng_st_save=leng_final;
leng_st_save2=leng_st_save./wei_final; % take the average if there are multiple streamlines btw two nodes
leng_st_save3=leng_st_save2*0.62; % convert to mm scale (0.62: average spacing btw tract segments)
leng_st_save3(isnan(leng_st_save3)==1)=0;

% weight matrix
wei_st_save=leng_final; % weight proportionally to the tract-length
wei_st_save2=wei_st_save./max(max(wei_st_save)); % normalize by dividing by maximum value (range: 0 - 1)

% ignore weight values less than 0.05 to reduce the storage size of the matrices
wei_st_save_11=wei_st_save2; wei_st_save_11(wei_st_save_11<0.05)=0; 
leng_st_save_11=leng_st_save3; leng_st_save_11(~wei_st_save_11)=0; 

save('patient folder path.../weights.mat','wei_st_save_11');
save('patient folder path.../tract_lengths.mat','leng_st_save_11');
