% cortical/subcortical surfaces reconstruction
% after running the code below, save the variables tri_s2, ver_s2, tri_c2, ver_c2 as .mat files in each patient folder

clear all

% load cortical/subcortical surfaces from TVB reconstruction output (surface_cort, surface_subcort)

ver_c=load('patient folder path.../surface_cort/vertices.txt');
ver_s=load('patient folder path.../surface_subcort/vertices.txt');
tri_c=load('patient folder path.../surface_cort/triangles.txt'); tri_c=tri_c+1; % python to matlab index
tri_s=load('patient folder path.../surface_subcort/triangles.txt'); tri_s=tri_s+1;
reg_map_c=load('patient folder path.../region_mapping_cort.destrieux.txt'); reg_map_c=reg_map_c+1; % python to matlab index
reg_map_s=load('patient folder path.../region_mapping_subcort.destrieux.txt'); reg_map_s=reg_map_s+1;

% ignore 0 index in region mapping files (excluding brainstem and some others) 

ind1_s=find(reg_map_s==0);
ver_s2=ver_s; ver_s2(ind1_s,:)=[];
tri_s2=tri_s;
for m1_s=1:1:length(ind1_s)
    ind2_s=ismember(tri_s2,ind1_s(m1_s));
    ind3_s=sum(ind2_s,2);
    ind4_s=find(ind3_s~=0);
    for m2_s=1:1:length(ind4_s)
        tri_s2(ind4_s(m2_s),:)=[-1 -1 -1];
    end
end    
ind5_s=ismember(tri_s2,[-1 -1 -1],'rows');
tri_s2(ind5_s,:)=[];

for m3_s=1:1:length(ver_s)
    ind6_s=length(find(ind1_s<m3_s));
    if ind6_s>0 && sum(sum(ismember(tri_s2,m3_s)))~=0
        ind7_s=ismember(tri_s2,m3_s);
        tri_s2(ind7_s)=tri_s2(ind7_s)-ind6_s;
    end
%     m3_s
end

ind1_c=find(reg_map_c==0);
ver_c2=ver_c; ver_c2(ind1_c,:)=[];
tri_c2=tri_c;
for m1_c=1:1:length(ind1_c)
    ind2_c=ismember(tri_c2,ind1_c(m1_c));
    ind3_c=sum(ind2_c,2);
    ind4_c=find(ind3_c~=0);
    for m2_c=1:1:length(ind4_c)
        tri_c2(ind4_c(m2_c),:)=[-1 -1 -1];
    end
end    
ind5_c=ismember(tri_c2,[-1 -1 -1],'rows');
tri_c2(ind5_c,:)=[];

for m3_c=1:1:length(ver_c)
    ind6_c=length(find(ind1_c<m3_c));
    if ind6_c>0 && sum(sum(ismember(tri_c2,m3_c)))~=0
        ind7_c=ismember(tri_c2,m3_c);
        tri_c2(ind7_c)=tri_c2(ind7_c)-ind6_c;
    end
    if mod(m3_c,100)==1
%         m3_c
    end
end

save('patient folder path.../tri_ig0_cort.mat','tri_c2');
save('patient folder path.../tri_ig0_subcort.mat','tri_s2');
save('patient folder path.../ver_ig0_cort.mat','ver_c2');
save('patient folder path.../ver_ig0_subcort.mat','ver_s2');
