function [ligand_new_str_structure,ligand_new_str_dG,ligand_new_str_repulse_num] = strEsearch_NRG_water(ligand_input,pocket_new)
tic
%ligand_input = lig_sel_modi;

%dG_select_passA = dG_select_pass_X1(index,1);
%repulse_num_loop = ligand_select_repulse_num_X1;



E_evaluate=0;
degree = -pi/6:pi/36:pi/6;
RcutoffPL = 6;
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;

const=halva*2;
%dG_select_pass=0;
ligand_inputA = ligand_input;



loop_num = 30;
Rmax = 0.4;
zero_coor = zeros(1,3);

clear dG_select

%if dG_select_pass-E_evaluate<=-0.5 & repulse_num_loopA>0

ln=1;
clear ligand_new_all
%cd('/home/dell/Octave_work/GARF_MT/GARF_MT/GARF_potential/MT_All_Contact/potential_hessian/');
ligand_output = cartvecNRG_waterF(pocket_new,ligand_inputA,loop_num);
%cd('/home/dell/Octave_work/Heatmap_Docking/');

ligand_new_str_structure = ligand_output.coordinate;
ligand_new_str_dG = ligand_output.dG;
ligand_new_str_repulse_num = ligand_output.repulse_num;
%ligand_new_str = ligand_input;

toc
