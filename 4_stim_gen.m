% select engaged fiber-tracts according to each active contact 
% and generate the stimulus input matrix accordingly

% enter the electrode contact positions in T1 space (in this study, 8 contacts, 4 each in the left and right hemispheres)

elec=[13.35 65.68 5.96; 14.35 65.80 9.31; 15.35 65.93 12.66; 16.34 66.05 16.02;...
      0.98 68.20 6.14; 0.63 68.33 9.62; 0.29 68.46 13.10; -0.06 68.60 16.58];

% load fiber-tracts streamlines from TVB reconstruction output
% (TVB recon pipeline reconstructs streamlines (coordinates) based on DTI space. Thus, these should be converted into T1 space using Flirt/?? commands in FSL, then loaded.)

for i1=0:1:99999
    file_n=sprintf('patient folder path.../SL_100k_T1_space/stream_lines-%d.txt',i1);
    data{i1+1}=dlmread(file_n,'',1,0);
end    
stream_lines=data;

% set the grid space according to electrode contact position 
% and load voltage distribution matrix (calculated from volt_dist_cal.m)

for k1=1:length(elec)
    for k2=1:3
        elec2=elec(k1,k2);
        if elec2-floor(elec2)>=0.75
            elec_d=1;
        elseif elec2-floor(elec2)>=0.25 && elec2-floor(elec2)<0.75
            elec_d=0.5;
        else
            elec_d=0;
        end
        elec_out(k1,k2)=floor(elec2)+elec_d;
    end
end

minX1 = min(elec_out(:,1))-10; maxX1 = max(elec_out(:,1))+10;
minY1 = min(elec_out(:,2))-10; maxY1 = max(elec_out(:,2))+10;
minZ1 = min(elec_out(:,3))-10; maxZ1 = max(elec_out(:,3))+10;
x = minX1:0.5:maxX1;
y = minY1:0.5:maxY1;
z = minZ1:0.5:maxZ1;
[xx, yy, zz] = meshgrid(x,y,z);
hx = x(2) - x(1); hy = y(2) - y(1); hz = z(2) - z(1);
Kx = (hy^2*hz^2)/(2*(hx^2*hy^2+hy^2*hz^2+hz^2*hx^2)); 
Ky = (hz^2*hx^2)/(2*(hx^2*hy^2+hy^2*hz^2+hz^2*hx^2));
Kz = (hx^2*hy^2)/(2*(hx^2*hy^2+hy^2*hz^2+hz^2*hx^2)); 

V_mat=load('patient folder path.../volt_dis_R0_contact.mat'); V=V_mat.V;

% derive how much voltage is applied to each segment for each fiber tract (streamline)
% and extract(select) the fiber-tracts if the voltage of the tract segment exceeds the certain threshold

m3=1; trk_act2=[]; 
for m1=1:length(stream_lines)
    trk=stream_lines{m1};
    trk2=ls_ch(trk);
    trk_act=[];
    for m2=1:size(trk2,1)
        trk_seg=trk2(m2,:);
        if ((trk_seg(1)>minX1)&&(trk_seg(1)<maxX1)) && ((trk_seg(2)>minY1)&&(trk_seg(2)<maxY1)) && ((trk_seg(3)>minZ1)&&(trk_seg(3)<maxZ1))
            a1=find(x==trk_seg(1));
            a2=find(y==trk_seg(2));
            a3=find(z==trk_seg(3));
        
            trk_act=[trk_act; [m2,V(a1,a2,a3)]];
        end
    end
    
    if isempty(trk_act)==0
        trk_act2{m3}=[m1*ones(size(trk_act,1),1),trk_act];
        m3=m3+1;
    end    
end

thre=1.5; % threshold voltage
n2=1; trk_act3=[];
for n1=1:length(trk_act2)
    trk_th=trk_act2{n1};
    trk_nz_ind=trk_th(1,1);
    if (sum(trk_th(:,3)>thre) > 0  || sum(trk_th(:,3)<-thre) > 0 ) 
        trk_act3{n2}=trk_th;
        n2=n2+1;
    end
end

% save the peak voltage of the selected fiber-tracts, their location, and the distances to both ends
% tr_act_info_th {}: [peak voltage, distances to both ends, coordinates of both ends of selected fiber-tracts]

str2=stream_lines;
tr_act_th=trk_act3;

tr_act_info_th=[];
for m1=1:length(tr_act_th)
    tr1=tr_act_th{m1}; tr1_ind=tr1(1,1);
    tr2=str2{tr1_ind}; tr2_len=length(tr2);
    tr3=[zeros(tr1(1,2)-1,1); tr1(:,3); zeros(tr2_len-tr1(end,2),1)];
    smooth_sig=-smoothdata(tr3,'gaussian',10);
    
    th_info=[];
    if smooth_sig(end)>smooth_sig(end-1)
        [pk,seg_loc]=findpeaks(smooth_sig);
        th_info=[pk,seg_loc];
        th_info=[th_info; [smooth_sig(end),length(smooth_sig)]];
    elseif smooth_sig(1)>smooth_sig(2)
        [pk,seg_loc]=findpeaks(smooth_sig);
        th_info=[pk,seg_loc];
        th_info=[th_info; [smooth_sig(1),1]];
    else
        [pk,seg_loc]=findpeaks(smooth_sig);
        th_info=[pk,seg_loc];
    end    
    
    ac_loc=th_info(:,2);    
    seg_act=find(smooth_sig(ac_loc)>=thre-0.1); 
    seg_act_val=smooth_sig(ac_loc(seg_act));
    seg_act_loc=ac_loc(seg_act);
           
    tr_info_br1=tr_act_th{m1};
    tr_info_br2=tr_info_br1(1,1);
    tr_info_pos1=str2{tr_info_br2};
    tr_info_dis_s2=[]; tr_info_dis_e2=[];
    
    for m2=1:length(seg_act_loc)
        tr_info_dis_s1=0.62*(seg_act_loc(m2)-1);
        tr_info_dis_e1=0.62*(size(tr_info_pos1,1)-(seg_act_loc(m2)));
        tr_info_dis_s2=[tr_info_dis_s2;tr_info_dis_s1];
        tr_info_dis_e2=[tr_info_dis_e2;tr_info_dis_e1];
    end
    
    tr_act_info_th{m1}=[seg_act_val,tr_info_dis_s2,tr_info_dis_e2];  
    if isempty(tr_act_info_th{m1})==1
        tr_act_info_th{m1}=tr_act_info_th{m1};
    else
        tr_pos_st_end=[tr_info_pos1(1,:),tr_info_pos1(end,:)];
        tr_act_info_th{m1}=[tr_act_info_th{m1},repelem(tr_pos_st_end,length(seg_act_loc),1)];
    end
end

% calculate the stimulus input to be applied to each node (based on selected fiber-tracts)

% (load cortical/subcortical surfaces from surf_const.m output)

tri_c2=load('patient folder path.../tri_ig0_cort.mat'); tri_c3=tri_c2.tri_c2;
tri_s2=load('patient folder path.../tri_ig0_subcort.mat'); tri_s3=tri_s2.tri_s2;
ver_c2=load('patient folder path.../ver_ig0_cort.mat'); ver_c3=ver_c2.ver_c2;
ver_s2=load('patient folder path.../ver_ig0_subcort.mat'); ver_s3=ver_s2.ver_s2;
[nf_c,nv_c] = reducepatch(tri_c3,ver_c3,1/30);
[nf_s,nv_s] = reducepatch(tri_s3,ver_s3,1/30);
ver_all=[nv_c;nv_s];
tri_all=[nf_c;nf_s+length(nv_c)];

speed_tr=6; 
stim_th=sparse(zeros(length(ver_all),500)); 
num_act_trk=0;
for m4=1:length(tr_act_info_th)
    tr_act1=tr_act_info_th{m4};
    if isempty(tr_act1)==0      
            for m5=1:size(tr_act1,1)
                tr_act2=tr_act1(m5,:);
                del_s=round(tr_act2(2)/speed_tr*5); % 1ms = 5 time steps 
                st_s=[zeros(1,del_s),ones(1,1)*((tr_act2(2)+tr_act2(3))/10)]; % weight proportionally to the tract-length
                st_s2=[st_s,zeros(1,size(stim_th,2)-length(st_s))];
                del_e=round(tr_act2(3)/speed_tr*5); 
                st_e=[zeros(1,del_e),ones(1,1)*((tr_act2(2)+tr_act2(3))/10)];
                st_e2=[st_e,zeros(1,size(stim_th,2)-length(st_e))];
                
                st_coor=tr_act2(1,4:6); end_coor=tr_act2(7:9);
                stim_th_st=sparse(zeros(length(ver_all),size(stim_th,2))); stim_th_end=sparse(zeros(length(ver_all),size(stim_th,2)));
                for v_coor=1:length(ver_all)
                    st_dis_v=sqrt(sum((ver_all(v_coor,:)-st_coor).^2));
                    if st_dis_v<10 
                        stim_th_st(v_coor,:)=stim_th_st(v_coor,:)+st_s2;
                    end                   
                    end_dis_v=sqrt(sum((ver_all(v_coor,:)-end_coor).^2));
                    if end_dis_v<10 
                        stim_th_end(v_coor,:)=stim_th_end(v_coor,:)+st_e2;
                    end 
                end
                
                if (length(nonzeros(stim_th_st)) > 0) && (length(nonzeros(stim_th_end)) > 0)                    
                    stim_th_each=stim_th_st+stim_th_end;
                    stim_th=stim_th+stim_th_each;     
                    num_act_trk=num_act_trk+1;
                end
            end   
    end
end

% save the generated stimulus input matrix 
dlmwrite('patient folder path.../stim_mat_R0.txt',full(stim_th),'delimiter','\t');

% visualize the selected fiber-tracts for the specific active contact
for m3=1:length(tr_act_th)
    tr1=tr_act_th{m3}; tr1_ind=tr1(1,1);
    tr_seg_coord=stream_lines{tr1_ind};
    
    figure(11); 
    plot3(tr_seg_coord(:,1),tr_seg_coord(:,2),tr_seg_coord(:,3),'k'); axis([-40 50 10 130 -30 60]); az =-90; el = 0; view(az, el); hold on;
    scatter3(elec(8,1),elec(8,2),elec(8,3),100,'filled');
end
hold off;
    