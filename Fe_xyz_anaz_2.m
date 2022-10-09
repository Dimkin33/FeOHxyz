clear,clc
close all
%load d10_50000.mat
load d16_30000.mat
%load d16_30000.mat
clearvars -except coord_xyz_f

%%

n_atom=19;

coord_Fe=coord_xyz_f(19,:);
coord_Fe(1:2)=0;
coord_xyz_f=coord_xyz_f-coord_Fe;
temp_zero=ones(1,19*30000*120)*0;
coord_xyz_f=[coord_xyz_f,temp_zero'];
c

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
    


%%
n_lmp=10;
time_n=30000;
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
    

