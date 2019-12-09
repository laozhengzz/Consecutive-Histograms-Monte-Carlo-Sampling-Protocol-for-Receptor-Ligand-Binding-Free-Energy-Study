function [ligand_final] = DioP_anchor_vf(ligand,Grid_Heatmap,pocket_centroid,protein,protein_centroidA)
%[ligand_final] =             DioP_anchor(ligand,Grid_Heatmap,pocket_centroid,pocket,protein);
%clear all
%load Grid_Heatmap_test_set
%tic
ligand_final=struct([]);
run GARF_Potential
RcutoffPL=6;
R_srch=0.5;
MTsample = R_srch;
R_L_tor_srch=0.1;
rd = 0.005;
hnstate = R_srch/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2;


ligand_centroid(1,1)=mean(ligand(:,1));
ligand_centroid(1,2)=mean(ligand(:,2));
ligand_centroid(1,3)=mean(ligand(:,3));


ligand_types = ligand(:,4);
ligand_types_all = ligand_types;
ligand_types = unique(ligand_types);


C_3=1;
C_2=2;
C_1=3;
C_ar=4;
O_3=5;
O_3p=6;
O_2=7;
O_co2=8;
O_2v=9;
N_2=10;
N_am=11;
N_pl3=12;
N_4=13;
P=14;
F=15;
Cl=16;
Br=17;
I=18;
C_cat=19;
S_3=20;
S_o=21;
HNa=22;
HOa=23;
HSa=24;
H3 =25;
C_3Oa=26;
C_3Na=27;
C_3La=28;
C_2Oa=29;
C_2Na=30;
C_2La=31;
C_arOa=32;
C_arNa=33;
C_arLa=34;
HCa=35;

priority_atomtype_list=[8	10
    7	10
    9	10
    13	10
    12	10
    5	10
    11	10
    10	10
    16	5
    17	5
    18	5
    6	2
    32	2
    33	2
    34	2
    29	2
    30	2
    31	2
    26	2
    27	2
    28	2
    20	0
    15	0
    4	0
    3	0
    2	0
    1	0
    19	0
    14	0
    21	0
    ];

top_prio_list = priority_atomtype_list( priority_atomtype_list(:,2)==10,:);
clear ligand_pin ligand_pin_type ligand_pin_distance
an=1;
for i = 1:size(ligand,1)
    for j =1:size(priority_atomtype_list( priority_atomtype_list(:,2)==max(priority_atomtype_list(:,2)),:),1)
        if ligand(i,4) == priority_atomtype_list(j,1)
            ligand_pin(an,1)=ligand(i,7);
            ligand_pin_type(an,1)=ligand(i,4);
            ligand_pin_distance(an,1)=sqrt((  ligand(ligand_pin(an,1),1) -ligand_centroid(1,1)        )^2+(     ligand(ligand_pin(an,1),2) -ligand_centroid(1,2)       )^2+(      ligand(ligand_pin(an,1),3) -ligand_centroid(1,3)       )^2);
            
            an=an+1;
            break
        end
    end
end
size(Grid_Heatmap,2)
size(ligand_pin_type,1)


clear E_grid_centroid_candi E_grid_centroid E_grid_boltzave
E_grid_boltzave = zeros(size(ligand_pin_type,1),size(Grid_Heatmap,2));
E_grid_centroid = zeros(size(Grid_Heatmap,2),3);
for i = 1:size(Grid_Heatmap,2)
    clear E_grid_centroid_candi RMSD_max
    E_grid_centroid_candi = zeros(size(ligand_pin_type,1),3);
    RMSD_max=zeros(size(ligand_pin_type,1),1);
    
    for j = 1:size(ligand_pin_type,1)
        eval( strcat  ('E_grid=Grid_Heatmap(i).Energy_',num2str(ligand_pin_type(j,1)),';'));
        if ~isempty(E_grid)
            %Egcount = floor(size(E_grid,1)/10);
            Egcount = size(E_grid,1);
            E_grid_centroid_candi(j,1) = mean(E_grid(1:Egcount,1));
            E_grid_centroid_candi(j,2) = mean(E_grid(1:Egcount,2));
            E_grid_centroid_candi(j,3) = mean(E_grid(1:Egcount,3));
            
            
            E_grid_boltzave(j,i) = boltzave(E_grid(:,4));
            
        end
    end
    
    E_grid_centroid_candi(E_grid_centroid_candi(:,1)==0 & E_grid_centroid_candi(:,2)==0 & E_grid_centroid_candi(:,3)==0,:)=[];
    E_grid_centroid(i,1)=mean(E_grid_centroid_candi(:,1));
    E_grid_centroid(i,2)=mean(E_grid_centroid_candi(:,2));
    E_grid_centroid(i,3)=mean(E_grid_centroid_candi(:,3));
    
end

dist_grid_cent = zeros(size(E_grid_centroid,1),size(pocket_centroid,1));

for i = 1:size(E_grid_centroid,1)
    if E_grid_centroid(i,1)~=0 && E_grid_centroid(i,2)~=0 && E_grid_centroid(i,3)~=0
        for j = 1:size(pocket_centroid,1)
            
            dist_grid_cent(i,j) = sqrt((  E_grid_centroid(i,1) -pocket_centroid(j,1)        )^2+(     E_grid_centroid(i,2) -pocket_centroid(j,2)       )^2+(      E_grid_centroid(i,3) -pocket_centroid(j,3)       )^2);
            
        end
    else
        dist_grid_cent(i,:) = 0;
    end
end

% dist_grid_cent: E_grid_centroid distances to the pocket centroids.

dist_grid_cent_nonzero = dist_grid_cent(dist_grid_cent(:,1)~=0,:);

Rmeasure=0;

clear ligand_select_X1 ligand_select_X2

if size(dist_grid_cent_nonzero,1)>=50
    Rmeasure=0.5;
elseif size(dist_grid_cent_nonzero,1)>=10 && size(dist_grid_cent_nonzero,1)<50
    Rmeasure=0.8;
elseif size(dist_grid_cent_nonzero,1)<10
    Rmeasure=1;
end

alm=1;
%clear ligand_final;
for pcnum = 1:size(pocket_centroid,1)
    %tic
    pocket_centroid_ind = pocket_centroid(pcnum,:);
    dist_pocket_cent = dist_grid_cent(:,pcnum);
    for lp =1:size(ligand_pin_distance,1)
        ligand_anchor1=ligand_centroid;
        ligand_anchor1(1,4:10)=0;
        ligand_anchor2=ligand(ligand_pin(lp,1),:);
        
        pockatomchosen_list=zeros(size(Grid_Heatmap,2),1);
        pac=1;
        for lg = 1:size(dist_pocket_cent,1)
            
            if abs(dist_pocket_cent(lg,1)-ligand_pin_distance(lp,1))<=Rmeasure
                pockatomchosen_list(pac,1)=lg;
                pac=pac+1;
            end
        end
        pockatomchosen_list(pockatomchosen_list(:,1)==0,:)=[];
        
        for j=1:size(pockatomchosen_list,1)
            %for j=2:2
            eval( strcat  ('grid_group1=Grid_Heatmap(',num2str(pockatomchosen_list(j,1)),').Energy_',num2str(ligand_pin_type(lp,1)),';'));
            
            grid_group1 = sortrows(grid_group1,[4]);
            grid_group1 = grid_group1(1:ceil(size(grid_group1,1)/3),:);
            grid_bin1 = pocket_centroid_ind;
            grid_bin2 = grid_filter3P(grid_group1);
            
            lgn=1;
            
            clear ligand_new_X1 ligand_new_X2 index_record dG_X
            dG_X_score=[];
            
            for gb1=1:size(grid_bin1,1)
                
                clear ligand_anchor1_S1 ligand_anchor2_S1 ligand_anchor3_S1 ligand_anchor1_S2
                ligand_anchor1_S1(1,1:3) = grid_bin1(gb1,1:3);
                ligand_anchor1_S1(1,4:10) = ligand_anchor1(1,4:10);
                
                ligand_anchor1_S2 = ligand_anchor1_S1;
                
                ligand_anchor2_S1(1,1) = ligand_anchor2(1,1)+(grid_bin1(gb1,1)-ligand_anchor1(1,1));
                ligand_anchor2_S1(1,2) = ligand_anchor2(1,2)+(grid_bin1(gb1,2)-ligand_anchor1(1,2));
                ligand_anchor2_S1(1,3) = ligand_anchor2(1,3)+(grid_bin1(gb1,3)-ligand_anchor1(1,3));
                ligand_anchor2_S1(1,4:10) = ligand_anchor2(1,4:10);
                
                ligand_S1(:,1) = ligand(:,1)+(grid_bin1(gb1,1)-ligand_anchor1(1,1));
                ligand_S1(:,2) = ligand(:,2)+(grid_bin1(gb1,2)-ligand_anchor1(1,2));
                ligand_S1(:,3) = ligand(:,3)+(grid_bin1(gb1,3)-ligand_anchor1(1,3));
                ligand_S1(:,4:10) = ligand(:,4:10);
                
                
                
                
                
                
                
                
                
                [ligand_S22,ligand_anchor2_S2] = molstand(ligand_S1,ligand_anchor1_S1,ligand_anchor2_S1,grid_bin1,gb1);
                
                
                for gb2 = 1:size(grid_bin2,1)
                %for gb2 = 13:13
                    clear ligand_anchor3_S2 ligand_anchor1_S3 ligand_anchor2_S3 ligand_anchor3_S3 ligand_anchor_S3 X1a X2a ligand_rot
                    
                    
                    
                    anglet = 0;
                    
                    for rotnum=1:24
                        
                        anglet = anglet+ pi/12;
                        [ligand_S2,ligand_anchor2_S22] = molrotatet(ligand_S22,ligand_anchor1_S2,ligand_anchor2_S2,anglet);
                        
                        %target = grid_bin2(gb2,:);
                        
                        %%%ligand_S3 = rotstrp2p(ligand_S2,ligand_anchor1_S2,ligand_anchor2_S22,target);
                        for axis = 1:3
                            [ligand_S2,ligand_anchor2_S22] = rotstrp2paxis(axis,grid_bin1(gb1,:),ligand_anchor2_S22,grid_bin2(gb2,:),ligand_S2);
                        end
                        
                        %ligand_S3 = molrotate(ligand_S2,ligand_anchor2_S22,grid_bin1,gb1,grid_bin2,gb2);
                        ligand_S3 = ligand_S2;
                        
                        ligand_anchor2_S3 = ligand_S3(ligand_S3(:,7)==ligand_anchor2_S22(1,7),:);
                        
                        dist_anchor2_gridbin2(rotnum,1) = dist(ligand_anchor2_S3,grid_bin2(gb2,:));
                        
                        
                        ligand_rot(rotnum).structure = ligand_S3;
                        
                        [protein_starterA_X1,protein_starterB_X1,repulse_num_X1] = pocket2find_PL_Dock(protein,ligand_rot(rotnum).structure,RcutoffPL);
                        [Part_Matr_com1_X1,CNN1_X1,VN1_X1] = inter_potential_PL_GARF_simple_A(protein_starterA_X1,protein_starterB_X1,vdw_d); %inter Protein-ligand interactions
                        
                        ZE_part_com1_X1 = -0.5918*Part_Matr_com1_X1 - VN1_X1*-0.5918;
                        ZN_part_com1_X1 = ((sqrt(CNN1_X1*8+1)+1)/2-4)*log(2)+((sqrt(CNN1_X1*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1_X1*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1_X1*8+1)+1)/2-1)*log(const);
                        dG_X1_sub(rotnum,1) = (ZE_part_com1_X1 + ZN_part_com1_X1*-0.5918)/20+3;
                        
                    end
                    
                    
                    maxrotind1 = find(dG_X1_sub == min(dG_X1_sub));
                    dG_X_score(lgn,1)=dG_X1_sub(maxrotind1(1,1),1);
                    dG_X(lgn).structure=ligand_rot(maxrotind1(1,1)).structure;
                    lgn=lgn+1;
                end
            end
            maxdgx=[];
            maxdgx = find(dG_X_score == min(dG_X_score));
            if ~isempty(maxdgx) && repulsenum(protein,dG_X(maxdgx(1,1)).structure)<2*size(ligand,1)
                ligand_final(alm).dG=dG_X_score(maxdgx(1,1),1);
                ligand_final(alm).structure=dG_X(maxdgx(1,1)).structure;
                ligand_final(alm).repulse=repulsenum(protein,dG_X(maxdgx(1,1)).structure);
                clear output_struct_centroid
                output_struct_centroid(1,1) = mean(dG_X(maxdgx(1,1)).structure(:,1));
                output_struct_centroid(1,2) = mean(dG_X(maxdgx(1,1)).structure(:,2));
                output_struct_centroid(1,3) = mean(dG_X(maxdgx(1,1)).structure(:,3));
                ligand_final(alm).procentroid_dist=dist(output_struct_centroid,protein_centroidA);
                alm=alm+1;
                
            end
        end
    end
    %toc
end












