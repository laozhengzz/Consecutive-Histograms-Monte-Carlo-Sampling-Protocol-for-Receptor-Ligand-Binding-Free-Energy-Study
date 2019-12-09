function rmsd_mark = RMSD_check(ligand_selected,ligand_candi,index_num)

rmsd_mark = 1;
for i = 1:index_num
    ligand = ligand_selected(index_num).structure;
    rmsd = RMSD_cal(ligand,ligand_candi);
    if rmsd<=2;
        rmsd_mark = 0;
        break
    end
end














