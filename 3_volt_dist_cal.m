% voltage distribution calculation according to active contact positions
% (based on finite difference method)
% after running the code below, save the variable V as .mat files in each patient folder

clear all

% enter the electrode contact positions in T1 space (in this study, 8 contacts, 4 each in the left and right hemispheres)

elec=[13.35 65.68 5.96; 14.35 65.80 9.31; 15.35 65.93 12.66; 16.34 66.05 16.02;...
      0.98 68.20 6.14; 0.63 68.33 9.62; 0.29 68.46 13.10; -0.06 68.60 16.58];

% round the contact positions to set the step size to 0.5 mm

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

% set the step size of mesh grids for FDM to 0.5mm
% calculate voltage distribution within the 10 mm boundary space based on the contact location

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
         
V = zeros(length(x),length(y),length(z));  

Vr1 = -8; Vr2 = 0; Vr3 = 0; Vr4 = 0; % apply voltage to the active contact 
Vl1 = 0; Vl2 = 0; Vl3 = 0; Vl4 = 0;
indxr1 = find(x==elec_out(1,1)); indyr1 = find(y==elec_out(1,2)); indzr1 = find(z==elec_out(1,3));
indxr2 = find(x==elec_out(2,1)); indyr2 = find(y==elec_out(2,2)); indzr2 = find(z==elec_out(2,3));
indxr3 = find(x==elec_out(3,1)); indyr3 = find(y==elec_out(3,2)); indzr3 = find(z==elec_out(3,3));
indxr4 = find(x==elec_out(4,1)); indyr4 = find(y==elec_out(4,2)); indzr4 = find(z==elec_out(4,3));
V(indxr1-1:indxr1+1,indyr1-1:indyr1+1,indzr1-1:indzr1+1) = Vr1; V(indxr2-1:indxr2+1,indyr2-1:indyr2+1,indzr2-1:indzr2+1) = Vr2;
V(indxr3-1:indxr3+1,indyr3-1:indyr3+1,indzr3-1:indzr3+1) = Vr3; V(indxr4-1:indxr4+1,indyr4-1:indyr4+1,indzr4-1:indzr4+1) = Vr4;
indxl1 = find(x==elec_out(5,1)); indyl1 = find(y==elec_out(5,2)); indzl1 = find(z==elec_out(5,3));
indxl2 = find(x==elec_out(6,1)); indyl2 = find(y==elec_out(6,2)); indzl2 = find(z==elec_out(6,3));
indxl3 = find(x==elec_out(7,1)); indyl3 = find(y==elec_out(7,2)); indzl3 = find(z==elec_out(7,3));
indxl4 = find(x==elec_out(8,1)); indyl4 = find(y==elec_out(8,2)); indzl4 = find(z==elec_out(8,3));
V(indxl1-1:indxl1+1,indyl1-1:indyl1+1,indzl1-1:indzl1+1) = Vl1; V(indxl2-1:indxl2+1,indyl2-1:indyl2+1,indzl2-1:indzl2+1) = Vl2;
V(indxl3-1:indxl3+1,indyl3-1:indyl3+1,indzl3-1:indzl3+1) = Vl3; V(indxl4-1:indxl4+1,indyl4-1:indyl4+1,indzl4-1:indzl4+1) = Vl4;

dSum = 1; n = 0; 

tol=0.001;

while  dSum > tol
    sum1 =  sum(sum(sum(V.^2)));
        
    for nx = 2: length(x)-1            
        for ny = 2: length(y)-1
            for nz = 2: length(z)-1                
                
                 V(indxr1-1:indxr1+1,indyr1-1:indyr1+1,indzr1-1:indzr1+1) = Vr1; % active contact
%                 V(indxr2-1:indxr2+1,indyr2-1:indyr2+1,indzr2-1:indzr2+1) = Vr2;
%                 V(indxr3-1:indxr3+1,indyr3-1:indyr3+1,indzr3-1:indzr3+1) = Vr3; 
%                 V(indxr4-1:indxr4+1,indyr4-1:indyr4+1,indzr4-1:indzr4+1) = Vr4;
%                 V(indxl1-1:indxl1+1,indyl1-1:indyl1+1,indzl1-1:indzl1+1) = Vl1; 
%                 V(indxl2-1:indxl2+1,indyl2-1:indyl2+1,indzl2-1:indzl2+1) = Vl2;
%                 V(indxl3-1:indxl3+1,indyl3-1:indyl3+1,indzl3-1:indzl3+1) = Vl3; 
%                 V(indxl4-1:indxl4+1,indyl4-1:indyl4+1,indzl4-1:indzl4+1) = Vl4;
                
                V(nx,ny,nz) = Kx * (V(nx+1,ny,nz) + V(nx-1,ny,nz)) + Ky * (V(nx,ny+1,nz) + V(nx,ny-1,nz)) + Kz * (V(nx,ny,nz+1) + V(nx,ny,nz-1));
            end        
        end
    end
    
   sum2 =  sum(sum(sum(V.^2)));
   dSum = abs(sum2 - sum1);
%    n = n+1
%    dSum
end

save('patient folder path.../volt_dis_R0_contact.mat','V');
