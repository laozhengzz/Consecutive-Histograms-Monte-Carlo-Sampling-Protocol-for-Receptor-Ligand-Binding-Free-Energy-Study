function [ligand_S2,ligand_anchor2_S2] = molstand(ligand_S1,ligand_anchor1_S1,ligand_anchor2_S1,grid_bin1,gb1)

distS12 = dist(ligand_anchor1_S1,ligand_anchor2_S1);
grid_bin2(1,1) = ligand_anchor1_S1(1,1);
grid_bin2(1,2) = ligand_anchor1_S1(1,2);
grid_bin2(1,3) = ligand_anchor1_S1(1,3)-distS12;
grid_bin2(1,4:10) = ligand_anchor2_S1(1,4:10);


for ls =1:size(ligand_S1,1)
    if ligand_S1(ls,7)~= ligand_anchor2_S1(1,7)
        [lsb1,lsb2] = find(ligand_S1(:,7)==ligand_anchor2_S1(1,7));
        ligand_anchor3 = ligand_S1(ls,:);
    elseif ligand_S1(ls,7)== ligand_anchor2_S1(1,7) && ls < size(ligand_S1,1)
        ligand_anchor3 = ligand_S1(ls+1,:);
        %continue
    elseif ligand_S1(ls,7)== ligand_anchor2_S1(1,7) && ls == size(ligand_S1,1)
        ligand_anchor3 = ligand_S1(ls-1,:);
        %continue
    end
    ligand_anchor3_S1 = ligand_anchor3;
    
    
    % rotate ligand_anchors_S1 according to grid_bin2
    
    
    rota=[ligand_anchor2_S1(1,1),ligand_anchor2_S1(1,2),ligand_anchor2_S1(1,3)];
    rotd=[grid_bin2(1,1),grid_bin2(1,2),grid_bin2(1,3)];
    rot0=[grid_bin1(gb1,1),grid_bin1(gb1,2),grid_bin1(gb1,3)];
    
    % x=1,y=2,z=3
    axis=3; % rotate along z axis
    [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rota,rotd,ligand_anchor3_S1);
    
    axis=1; % rotate along y axis
    
    [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rotaN,rotd,ligand_anchor3_rotaN);
    
    axis=2; % rotate along x axis
    [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rotaN,rotd,ligand_anchor3_rotaN);
    
    ligand_anchor2_S2(1,1:3)=rotaN(1,1:3);
    ligand_anchor2_S2(1,4:10)=ligand_anchor2_S1(1,4:10);
    ligand_anchor3_S2(1,1:3)=ligand_anchor3_rotaN(1,1:3);
    ligand_anchor3_S2(1,4:10)=ligand_anchor3(1,4:10);
    
    if ligand_S1(ls,7)~= ligand_anchor2_S1(1,7)
        ligand_S2(lsb1,:)=ligand_anchor2_S2;
        ligand_S2(ls,:)=ligand_anchor3_S2;
    elseif ligand_S1(ls,7)== ligand_anchor2_S1(1,7) && ls < size(ligand_S1,1)
        ligand_S2(ls,:)=ligand_anchor2_S2;
        ligand_S2(ls+1,:)=ligand_anchor3_S2;
    elseif ligand_S1(ls,7)== ligand_anchor2_S1(1,7) && ls == size(ligand_S1,1)
        ligand_S2(ls,:)=ligand_anchor2_S2;
        ligand_S2(ls-1,:)=ligand_anchor3_S2;
        
        
    end
    
end