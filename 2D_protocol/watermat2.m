function water_mat = watermat2(pocket_structure,pocket_centroidA,pocket_centroid,lcp_dist_vec)
%pocket_structure = protein;
op_radi=[0
    3.325
3.73
3.65
3.635
3.66
3.155
2.895
2.85
2.955
2.785
2.84
2.845
2.735
2.82
2.745
3.015
2.76
3.78
];
grid_d=3.2;
max_x=0+10;
min_x=0-10;
max_y=0+10;
min_y=0-10;
max_z=0+10;
min_z=0-10;
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

P_sh_grid_stand=grid_coor;

clear P_sh_grid
%E_grid=zeros(10^7,1);
P_sh_grid(:,1) = P_sh_grid_stand(:,1) + pocket_centroidA(1,1);
P_sh_grid(:,2) = P_sh_grid_stand(:,2) + pocket_centroidA(1,2);
P_sh_grid(:,3) = P_sh_grid_stand(:,3) + pocket_centroidA(1,3);




for jz=1:1:size(pocket_structure,1)
%for jz=1:1:1
    %tic
    %for jz=1:1:1
    %if pocket_structure(jz,4)~=8
        %clear P_sh_grid
        
        %P_sh_grid(:,1) = P_sh_grid_stand(:,1) + pocket_structure(jz,1);
        %P_sh_grid(:,2) = P_sh_grid_stand(:,2) + pocket_structure(jz,2);
        %P_sh_grid(:,3) = P_sh_grid_stand(:,3) + pocket_structure(jz,3);
        
        
        R_point=pocket_structure(jz,:);
        lig_blackwhole=pocket_structure(sqrt( (pocket_structure(:,1)-R_point(1,1)).^2+(pocket_structure(:,2)-R_point(1,2)).^2+(pocket_structure(:,3)-R_point(1,3)).^2 )<5.6,:);
        
        size_lig_blackwhole=size(lig_blackwhole);
        for jj=1:size_lig_blackwhole(1,1)
            if lig_blackwhole(jj,7)~=pocket_structure(jz,7)
                P_sh_grid(sqrt( (P_sh_grid(:,1)-lig_blackwhole(jj,1)).^2+(P_sh_grid(:,2)-lig_blackwhole(jj,2)).^2+(P_sh_grid(:,3)-lig_blackwhole(jj,3)).^2 )<=lig_blackwhole(jj,6),:)=[];
            end
        end
        
        
    %end
    %toc
end

water_mat=[];
for i = 1:size(lcp_dist_vec,1)
    water_mat_temp = P_sh_grid;
    
    water_mat_temp(sqrt( (water_mat_temp(:,1)-pocket_centroid(1,1)).^2+(water_mat_temp(:,2)-pocket_centroid(1,2)).^2+(water_mat_temp(:,3)-pocket_centroid(1,3)).^2 )<lcp_dist_vec(i,1)-0.5,:)=[];
    water_mat_temp(sqrt( (water_mat_temp(:,1)-pocket_centroid(1,1)).^2+(water_mat_temp(:,2)-pocket_centroid(1,2)).^2+(water_mat_temp(:,3)-pocket_centroid(1,3)).^2 )>lcp_dist_vec(i,1)+3.0,:)=[];
    water_mat(size(water_mat,1)+1:size(water_mat,1)+size(water_mat_temp,1),:) = water_mat_temp;
end

water_mat = unique(water_mat,'rows');

water_mat(:,4) = 20;
water_mat(:,5) = 16;
water_mat(:,6) = op_radi(16,1);

%water_num=[1:1:size(water_mat,1)]';

%water_mat(:,7)=water_num;
water_mat(:,8)=1;
water_mat(:,9)=1;
water_mat(:,10)=101;
















