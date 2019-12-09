function [protein_surface] = pocketsurface(lyr,centroid,protein_centroidA,pocket_protein_centroid_dist,protein_structure,grid_d)
   
   %centroid = pocket_centroidA;
   %protein_structure = protein;
   
   P_sh_grid=[];
   
%grid_d=0.4;
Rmax = lyr+2.61;
Rmax2 = lyr+2*grid_d;
max_x=0+Rmax2;
min_x=0-Rmax2;
max_y=0+Rmax2;
min_y=0-Rmax2;
max_z=0+Rmax2;
min_z=0-Rmax2;
sig_x=transpose(min_x:grid_d:max_x);
sig_y=transpose(min_y:grid_d:max_y);
sig_z=transpose(min_z:grid_d:max_z);

size_sig_x=size(sig_x);
size_sig_y=size(sig_y);
size_sig_z=size(sig_z);


for q=1:1:size_sig_x(1,1)
    grid_coor((q-1)*size_sig_y(1,1)*size_sig_z(1,1)+1:q*size_sig_y(1,1)*size_sig_z(1,1),1)=sig_x(q,1);
end
%toc
%tic
for q=1:1:size_sig_y(1,1)
    grid_coor1r((q-1)*size_sig_z(1,1)+1:q*size_sig_z(1,1),1)=sig_y(q,1);
end
%toc
%tic
grid_coor(:,2)=repmat(grid_coor1r,size_sig_x(1,1),1);
grid_coor(:,3)=repmat(sig_z,size_sig_x(1,1)*size_sig_y(1,1),1);
P_sh_num=0;
P_sh_grid_stand = zeros(10^5,3);
%for lyr = 2.7:grid_d:3.8
%for lyr = 2.7:grid_d:4    
    %clear P_sh
    P_sh=grid_coor(  (sqrt( (grid_coor(:,1)-0).^2+(grid_coor(:,2)-0).^2+(grid_coor(:,3)-0).^2 )<= lyr+grid_d/2)&(sqrt( (grid_coor(:,1)-0).^2+(grid_coor(:,2)-0).^2+(grid_coor(:,3)-0).^2 )>= lyr-grid_d/2),:);
    P_sh_grid_stand(P_sh_num+1:P_sh_num+size(P_sh,1),1:3)=P_sh;
    P_sh_num = P_sh_num+size(P_sh,1);
%end

P_sh_grid_stand(P_sh_grid_stand(:,1)==0,:)=[];
P_sh_grid_stand = unique(P_sh_grid_stand,'rows');

P_sh_grid(:,1) = P_sh_grid_stand(:,1) + centroid(1,1);
P_sh_grid(:,2) = P_sh_grid_stand(:,2) + centroid(1,2);
P_sh_grid(:,3) = P_sh_grid_stand(:,3) + centroid(1,3);

protein_structure_temp = protein_structure;
protein_structure_temp( sqrt( (protein_structure_temp(:,1)-centroid(1,1)).^2+(protein_structure_temp(:,2)-centroid(1,2)).^2+(protein_structure_temp(:,3)-centroid(1,3)).^2 )>Rmax,:)=[];


if size(protein_structure_temp,1)>0
    for jz=1:1:size(protein_structure_temp,1)
        P_sh_grid(sqrt( (P_sh_grid(:,1)-protein_structure_temp(jz,1)).^2+(P_sh_grid(:,2)-protein_structure_temp(jz,2)).^2+(P_sh_grid(:,3)-protein_structure_temp(jz,3)).^2 )< protein_structure_temp(jz,6),:)=[];
    end
end

max_pp_d = min(pocket_protein_centroid_dist);

P_sh_grid(sqrt( (P_sh_grid(:,1)-protein_centroidA(1,1)).^2+(P_sh_grid(:,2)-protein_centroidA(1,2)).^2+(P_sh_grid(:,3)-protein_centroidA(1,3)).^2 )< max_pp_d,:)=[];

%P_sh_grid_mark = zeros(size(P_sh_grid,1),1);
%for i = 1:size(P_sh_grid,1)
%    dist_vec = dist(P_sh_grid,P_sh_grid(i,:));
%    if ~isempty(dist_vec)
%        dist_vec(dist_vec(:,1)>1.6 | dist_vec(:,1)==0,:)=[];
        
%        if size(dist_vec,1)<=10
%            P_sh_grid_mark(i,1)=1;
%        end
%    end
    
%end
    
%P_sh_grid(P_sh_grid_mark(:,1)==1,:)=[];



%P_sh_grid(sqrt( (P_sh_grid(:,1)-P_sh_grid_ave(1,1)).^2+(P_sh_grid(:,2)-P_sh_grid_ave(1,2)).^2+(P_sh_grid(:,3)-P_sh_grid_ave(1,3)).^2 )> 5,:)=[];
protein_surface=[];
nps = 1;
for i = 1:size(protein_structure,1)
    dist_pg = dist(P_sh_grid,protein_structure(i,:));
    if min(dist_pg)<protein_structure(i,6)+grid_d;
        protein_surface(nps,:) = protein_structure(i,:);
        nps = nps+1;
    end
end
    
%protein_surface(protein_surface(:,4)==8,:)=[];

