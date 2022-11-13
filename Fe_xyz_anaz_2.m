clear,clc
close all
%load d10_50000.mat
%load d16_30000.mat
load ../datafile/d16_45000.mat
clearvars -except coord_xyz_f
time_n=45000;
n_lmp=120;

%%

n_atom=19;

coord_Fe=coord_xyz_f(19,:);
coord_Fe(1:2)=0;
coord_xyz_f=coord_xyz_f-coord_Fe;
temp_zero=ones(1,n_atom*time_n*n_lmp)*0;
coord_xyz_f=[coord_xyz_f,temp_zero'];


%% pair definition
clearvars temp_zero
coord_one=coord_xyz_f(1:n_atom,3:5);
temp_zero=ones(1,n_atom)*0;
coord_one=[coord_one,temp_zero'];
for i=7:n_atom-1
    [m,ind]=min(vecnorm(coord_one(i,1:3)-coord_one(1:6,1:3),2,2));
    [m,ind,i];
    coord_one(i,4)=ind;
end
    
%t2df

%%

for i=1:n_lmp*time_n
    for i_atom=1:6
        coord_xyz_f(i_atom+(i-1)*n_atom,6)=vecnorm(coord_xyz_f(i_atom+(i-1)*n_atom,3:5));
    end
    shift1=6;
    for i_atom=7:2:17
        coord_xyz_f(i_atom+(i-1)*n_atom,6)=vecnorm(coord_xyz_f(i_atom+(i-1)*n_atom,3:5)-coord_xyz_f(i_atom-shift1+(i-1)*n_atom,3:5));
        coord_xyz_f(i_atom+(i-1)*n_atom+1,6)=vecnorm(coord_xyz_f(i_atom+(i-1)*n_atom+1,3:5)-coord_xyz_f(i_atom-shift1+(i-1)*n_atom,3:5));
        shift1=shift1+1;
    end
end
  
%% angle section
for i=1:n_lmp*time_n
    %shiftl=6;
    for i_atom=1:6
        
        OFe=coord_xyz_f(i_atom+(i-1)*n_atom,3:5)-coord_xyz_f(19+(i-1)*n_atom,3:5);
        OH1=coord_xyz_f(i_atom+(i-1)*n_atom,3:5)-coord_xyz_f(2*i_atom+5+(i-1)*n_atom,3:5);
        OH2=coord_xyz_f(i_atom+(i-1)*n_atom,3:5)-coord_xyz_f(2*i_atom+6+(i-1)*n_atom,3:5);
        
        O2Fe=coord_xyz_f(mod(i_atom,6)+1+(i-1)*n_atom,3:5)-coord_xyz_f(19+(i-1)*n_atom,3:5);
        OOFe_normal=cross(OFe,O2Fe);
        OH_normal=cross(OH1,OH2);
        coord_xyz_f(i_atom+(i-1)*n_atom,7)=acos(dot(OFe,OH_normal)/(norm(OFe)*norm(OH_normal)))*180/pi;
        coord_xyz_f(i_atom+(i-1)*n_atom,8)=acos(dot(OH1,OOFe_normal)/(norm(OH1)*norm(OOFe_normal)))*180/pi;
        %shiftl=shiftl+1;
    end
end
%%
clearvars len_O len_O_0 O_time1
len_O=reshape(coord_xyz_f(:,6),n_atom,[]);
len_O(7:19,:)=[];
len_O=len_O';
O_time=mean(len_O,2);
O_time=reshape(O_time,[],n_lmp);
clearvars len_O


%% angle  norm(HOH) and OFe
len_angle=reshape(coord_xyz_f(:,7),n_atom,[]);
len_angle(7:19,:)=[];
len_angle=len_angle';
angle_time=mean(len_angle,2);
angle_time=reshape(angle_time,[],n_lmp);
clearvars len_angle

%% angle norm(OFeO) and OH
len_angle2=reshape(coord_xyz_f(:,8),n_atom,[]);
len_angle2(7:19,:)=[];
len_angle2=len_angle2';
angle_time2=mean(len_angle2,2);
angle_time2=reshape(angle_time2,[],n_lmp);
clearvars len_angle2


%%
time_plot=1000;
n_lmp_plot=[1, 30,120];
plot(1:time_plot,mean(O_time(1:time_plot,1:n_lmp_plot(1)),2),...
    1:time_plot,mean(O_time(1:time_plot,1:n_lmp_plot(2)),2),...
    1:time_plot,mean(O_time(1:time_plot,1:n_lmp_plot(3)),2));
     title('\rho, distance O (water) - Fe ion')
     xlabel('time, fs ')
     ylabel('\rho, Angtsrem')
     legend(strtrim(cellstr(num2str(n_lmp_plot'))'))
%%
figure
     time_plot=1000;
n_lmp_plot=[1, 30,120];
plot(1:time_plot,mean(angle_time(1:time_plot,1:n_lmp_plot(1)),2),...
    1:time_plot,mean(angle_time(1:time_plot,1:n_lmp_plot(2)),2),...
    1:time_plot,mean(angle_time(1:time_plot,1:n_lmp_plot(3)),2));
     title('\alpha, angle between norm(HOH) and OFe')
     xlabel('time, fs ')
     ylabel('\alpha, degree')
     legend(strtrim(cellstr(num2str(n_lmp_plot'))'))
%%
     figure
     time_plot=1000;
n_lmp_plot=[1, 30,120];
plot(1:time_plot,mean(angle_time2(1:time_plot,1:n_lmp_plot(1)),2),...
    1:time_plot,mean(angle_time2(1:time_plot,1:n_lmp_plot(2)),2),...
    1:time_plot,mean(angle_time2(1:time_plot,1:n_lmp_plot(3)),2));
     title('\beta, beta between norm(OFeO) and OH')
     xlabel('time, fs ')
     ylabel('\beta, degree')
     legend(strtrim(cellstr(num2str(n_lmp_plot'))'))
     
