function [P_grid_final,P_grid_final_ave] = gridshell_xy(lyr1,lyr2,Panchor1,Panchor2,ligand_centroid,protein_centroidA,protein_structure,pocket,grid_d,max_dist_llc)

%lyr1 = lyrx;
%lyr2 = lyrz;
%Panchor1 = Panchorx;
%Panchor2 = Panchorz;
%protein_structure = protein;

dpp = dist(ligand_centroid,protein_centroidA);


P_grid_final = [];
P_grid_final_ave=[];
%centroid = pocket_centroidA;
%protein_structure = protein;

P_sh_grid_ave = [];

P_sh_grid1=[];
P_sh_grid2=[];


P_sh_grid1 = gridlayer_axis(lyr1,grid_d,Panchor1,ligand_centroid,protein_structure,1);
P_sh_grid2 = gridlayer_axis(lyr2,grid_d,Panchor2,ligand_centroid,protein_structure,3);

%P_sh_grid1 = gridlayer(lyr1,grid_d,Panchor1,protein_structure);
%P_sh_grid2 = gridlayer(lyr2,grid_d,Panchor2,protein_structure);

%max_pp_d = min(pocket_protein_centroid_dist);

%P_sh_grid1(sqrt( (P_sh_grid1(:,1)-protein_centroidA(1,1)).^2+(P_sh_grid1(:,2)-protein_centroidA(1,2)).^2+(P_sh_grid1(:,3)-protein_centroidA(1,3)).^2 )< max_pp_d,:)=[];
%P_sh_grid2(sqrt( (P_sh_grid2(:,1)-protein_centroidA(1,1)).^2+(P_sh_grid2(:,2)-protein_centroidA(1,2)).^2+(P_sh_grid2(:,3)-protein_centroidA(1,3)).^2 )< max_pp_d,:)=[];


DistPPgrid=[];



for i = 1:size(P_sh_grid1,1)
    clear DistPPgrid_temp
    DistPPgrid_temp = dist(P_sh_grid2,P_sh_grid1(i,:));
    DistPPgrid_temp(:,2) = i;
    DistPPgrid_temp(:,3) = 1:size(DistPPgrid_temp,1);
    
    DistPPgrid_temp(DistPPgrid_temp(:,1) > grid_d/2,:)=[];
    
    DistPPgrid(size(DistPPgrid,1)+1:size(DistPPgrid,1)+size(DistPPgrid_temp,1),:) = DistPPgrid_temp;
end

for i = 1:size(DistPPgrid,1)
    P_grid_final(i,1) = (P_sh_grid1(DistPPgrid(i,2),1) + P_sh_grid2(DistPPgrid(i,3),1))/2 ;
    P_grid_final(i,2) = (P_sh_grid1(DistPPgrid(i,2),2) + P_sh_grid2(DistPPgrid(i,3),2))/2 ;
    P_grid_final(i,3) = (P_sh_grid1(DistPPgrid(i,2),3) + P_sh_grid2(DistPPgrid(i,3),3))/2 ;
end

mark_distgP = zeros(size(P_grid_final,1),1);
for i = 1:size(P_grid_final,1)
    clear distgP
    min_distgP = min(dist(pocket,P_grid_final(i,:)));
    if min_distgP > max_dist_llc + 2.5
        mark_distgP(i,1) = 1;
    end
end
    
P_grid_final(mark_distgP(:,1) == 1,:)=[];

mark_distgP = zeros(size(P_grid_final,1),1);
for i = 1:size(P_grid_final,1)
    clear distlcg distpcg
    
    distlcg = dist(ligand_centroid,P_grid_final(i,:));
    distpcg = dist(protein_centroidA,P_grid_final(i,:));
    if distlcg^2 > dpp^2 + distpcg^2
    %if distlcg^2 + dpp^2 > distpcg^2
        mark_distgP(i,1) = 1;
    end
end
    
P_grid_final(mark_distgP(:,1) == 1,:)=[];


if ~isempty(P_grid_final)
    P_grid_final_ave(1,1) = mean(P_grid_final(:,1));
    P_grid_final_ave(1,2) = mean(P_grid_final(:,2));
    P_grid_final_ave(1,3) = mean(P_grid_final(:,3));
end
%P_sh_grid(sqrt( (P_sh_grid(:,1)-P_sh_grid_ave(1,1)).^2+(P_sh_grid(:,2)-P_sh_grid_ave(1,2)).^2+(P_sh_grid(:,3)-P_sh_grid_ave(1,3)).^2 )> 5,:)=[];




