function [ligand_final_str,ligand_final_dG,ligand_final_repulse_num] = Erelax_NRG_water(ligand_select_X1,dG_select_pass_X1,index,protein,water_mat,vdw_d)

%index = i;

protein_score = [protein
    water_mat];

rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;
const=halva*2;
%index = i
ligand_select_coor_X1 = ligand_select_X1(dG_select_pass_X1(index,3)).coordinate;
ligand_select_cent_X1(1,1) = mean(ligand_select_coor_X1(:,1));
ligand_select_cent_X1(1,2) = mean(ligand_select_coor_X1(:,2));
ligand_select_cent_X1(1,3) = mean(ligand_select_coor_X1(:,3));

lig_sel_modi(:,1) = ligand_select_coor_X1(:,1)+ 0-ligand_select_cent_X1(1,1);
lig_sel_modi(:,2) = ligand_select_coor_X1(:,2)+ 0-ligand_select_cent_X1(1,2);
lig_sel_modi(:,3) = ligand_select_coor_X1(:,3)+ 0-ligand_select_cent_X1(1,3);
lig_sel_modi(:,4:10) = ligand_select_coor_X1(:,4:10);

pocket_new(:,1) = protein_score(:,1)+ 0-ligand_select_cent_X1(1,1);
pocket_new(:,2) = protein_score(:,2)+ 0-ligand_select_cent_X1(1,2);
pocket_new(:,3) = protein_score(:,3)+ 0-ligand_select_cent_X1(1,3);
pocket_new(:,4:10) = protein_score(:,4:10);

for i = 1:size(lig_sel_modi,1)
    pocket_new(pocket_new(:,10) == 101 & min(dist(pocket_new,lig_sel_modi(i,:)))<2,:)=[];
end




[ligand_final_str,ligand_relax_dG,ligand_final_repulse_num] = strEsearch_NRG_water(lig_sel_modi,pocket_new);
ligand_final_str(:,1) = ligand_final_str(:,1) + ligand_select_cent_X1(1,1);
ligand_final_str(:,2) = ligand_final_str(:,2) + ligand_select_cent_X1(1,2);
ligand_final_str(:,3) = ligand_final_str(:,3) + ligand_select_cent_X1(1,3);
ligand_final_str(:,4:10)=ligand_final_str(:,4:10);
%rmsd_new_strA = RMSD_cal(ligand,ligand_final_str)
RcutoffPL = 6;
[protein_starterA,protein_starterB,repulse_num] = pocket2find_PL_Dock_water(protein_score,ligand_final_str,RcutoffPL);
[Part_Matr_com1_new,CNN1_new,VN1_new] = inter_potential_PL_GARF_simple_water(protein_starterA,protein_starterB,vdw_d); %inter Protein-ligand interactions

ZE_part_com1_new = -0.5918*Part_Matr_com1_new - VN1_new*-0.5918;
ZN_part_com1_new = ((sqrt(CNN1_new*8+1)+1)/2-4)*log(2)+((sqrt(CNN1_new*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1_new*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1_new*8+1)+1)/2-1)*log(const);
ligand_final_dG = (ZE_part_com1_new + ZN_part_com1_new*-0.5918)/20+3;

%input_struc = ligand_new_strA;
%list = '1lri';
%output_dir = 'F:\1lri docking file\output';
%asdf = 2;
%out = lig_output(input_struc,list,asdf,output_dir);



%end


