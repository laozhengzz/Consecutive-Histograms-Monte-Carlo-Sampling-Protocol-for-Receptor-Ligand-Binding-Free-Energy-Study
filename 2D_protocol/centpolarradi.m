function dist_vec = centpolarradi(ligand_centroid,ligand)

ip=1;
for i = 1:size(ligand,1)
    if ligand(i,4) >= 5 && ligand(i,4) <= 18
        dist_vec(ip,1)=dist(ligand(i,:),ligand_centroid);
        ip=ip+1;
    end
end