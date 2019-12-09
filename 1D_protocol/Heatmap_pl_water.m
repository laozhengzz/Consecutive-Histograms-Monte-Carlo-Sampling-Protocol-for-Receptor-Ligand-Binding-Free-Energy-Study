function Grid_Heatmap = Heatmap_pl_water(pocket_structure,protein_whole,ligand)

ligand_types = ligand(:,4);
ligand_types = unique(ligand_types);

grid_d=0.2;
max_x=0+5;
min_x=0-5;
max_y=0+5;
min_y=0-5;
max_z=0+5;
min_z=0-5;
sig_x=transpose(min_x:grid_d:max_x);
sig_y=transpose(min_y:grid_d:max_y);
sig_z=transpose(min_z:grid_d:max_z);

size_sig_x=size(sig_x);
size_sig_y=size(sig_y);
size_sig_z=size(sig_z);


for q=1:1:size_sig_x(1,1)
    grid_coor((q-1)*size_sig_y(1,1)*size_sig_z(1,1)+1:q*size_sig_y(1,1)*size_sig_z(1,1),1)=sig_x(q,1);
end
%%toc
%%tic
for q=1:1:size_sig_y(1,1)
    grid_coor1r((q-1)*size_sig_z(1,1)+1:q*size_sig_z(1,1),1)=sig_y(q,1);
end
%%toc
%%tic
grid_coor(:,2)=repmat(grid_coor1r,size_sig_x(1,1),1);
grid_coor(:,3)=repmat(sig_z,size_sig_x(1,1)*size_sig_y(1,1),1);
P_sh_num=0;
P_sh_grid_stand = zeros(10^5,3);
%for lyr = 2.7:grid_d:3.8
for lyr = 2.7:grid_d:3.2    
    %clear P_sh
    P_sh=grid_coor(  (sqrt( (grid_coor(:,1)-0).^2+(grid_coor(:,2)-0).^2+(grid_coor(:,3)-0).^2 )<= lyr+grid_d)&(sqrt( (grid_coor(:,1)-0).^2+(grid_coor(:,2)-0).^2+(grid_coor(:,3)-0).^2 )> lyr),:);
    P_sh_grid_stand(P_sh_num+1:P_sh_num+size(P_sh,1),1:3)=P_sh;
    P_sh_num = P_sh_num+size(P_sh,1);
end

P_sh_grid_stand(P_sh_grid_stand(:,1)==0,:)=[];
P_sh_grid_stand = unique(P_sh_grid_stand,'rows');


pocket_structure(:,6) = pocket_structure(:,6)-1.6;
protein_whole(:,6) = protein_whole(:,6)-1.6;
run GARF_Potential_Para
for jz=1:1:size(pocket_structure,1)
%for jz=1:1:1
    %tic
    %for jz=1:1:1
    if pocket_structure(jz,4)~=8
        %%tic
        clear P_sh_grid
        %E_grid=zeros(10^7,1);
        P_sh_grid(:,1) = P_sh_grid_stand(:,1) + pocket_structure(jz,1);
        P_sh_grid(:,2) = P_sh_grid_stand(:,2) + pocket_structure(jz,2);
        P_sh_grid(:,3) = P_sh_grid_stand(:,3) + pocket_structure(jz,3);
        %grid_coor(:,4)=zz;
        %size_grid_coor=size(grid_coor);
        %%%%%%% grid_coor raw grid point box %%%%%%%%%
        
        R_point=pocket_structure(jz,:);
        lig_blackwhole=protein_whole(sqrt( (protein_whole(:,1)-R_point(1,1)).^2+(protein_whole(:,2)-R_point(1,2)).^2+(protein_whole(:,3)-R_point(1,3)).^2 )<5.6,:);
        
        size_lig_blackwhole=size(lig_blackwhole);
        for jj=1:size_lig_blackwhole(1,1)
            if lig_blackwhole(jj,7)~=pocket_structure(jz,7)
                %P_sh_grid(sqrt( (P_sh_grid(:,1)-lig_blackwhole(jj,1)).^2+(P_sh_grid(:,2)-lig_blackwhole(jj,2)).^2+(P_sh_grid(:,3)-lig_blackwhole(jj,3)).^2 )<=lig_blackwhole(jj,6),:)=[];
                P_sh_grid(sqrt( (P_sh_grid(:,1)-lig_blackwhole(jj,1)).^2+(P_sh_grid(:,2)-lig_blackwhole(jj,2)).^2+(P_sh_grid(:,3)-lig_blackwhole(jj,3)).^2 )<=2.6,:)=[];
            end
        end
        
        
        rd = 0.005;
        hnstate = grid_d/rd;
        halva=hnstate;
        const=halva*2+1;
        Dist_Mat_prep = zeros(size(P_sh_grid,1),const);
        Dist_Mat_prep(:,halva+1)=0;
        
        for j = 1:halva
            Dist_Mat_prep(:,halva+1-j)=0-j*rd;
            Dist_Mat_prep(:,halva+1+j)=j*rd;
        end
        
        for i =1:size(ligand_types,1)
            grid_type = ligand_types(i,1);
            E_grid = mtgrid(pocket_structure,protein_whole,P_sh_grid,grid_type,Paraset,grid_d,Dist_Mat_prep,jz);
            eval( strcat  ('Grid_Heatmap(jz).Energy_',num2str(grid_type),'=E_grid;'));
            
        end
        
    end
    %toc
end