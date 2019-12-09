function ligand_S3 = molrotate(ligand_S2,ligand_anchor2_S22,grid_bin1,gb1,grid_bin2,gb2)





for ls =1:size(ligand_S2,1)
    if ligand_S2(ls,7)~= ligand_anchor2_S22(1,7)
        [lsb1,lsb2] = find(ligand_S2(:,7)==ligand_anchor2_S22(1,7));
        ligand_anchor3 = ligand_S2(ls,:);
    
    ligand_anchor3_S1 = ligand_anchor3;
    
    
    % rotate ligand_anchors_S1 according to grid_bin2
    
    
    rota=[ligand_anchor2_S22(1,1),ligand_anchor2_S22(1,2),ligand_anchor2_S22(1,3)];
    rotd=[grid_bin2(gb2,1),grid_bin2(gb2,2),grid_bin2(gb2,3)];
    rot0=[grid_bin1(gb1,1),grid_bin1(gb1,2),grid_bin1(gb1,3)];
    
    % x=1,y=2,z=3
    axis=3; % rotate along z axis
    [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rota,rotd,ligand_anchor3_S1);
    
    axis=1; % rotate along y axis
    
    [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rotaN,rotd,ligand_anchor3_rotaN);
    
    axis=2; % rotate along x axis
    [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rotaN,rotd,ligand_anchor3_rotaN);
    
    ligand_anchor2_S2(1,1:3)=rotaN(1,1:3);
    ligand_anchor2_S2(1,4:10)=ligand_anchor2_S22(1,4:10);
    ligand_anchor3_S2(1,1:3)=ligand_anchor3_rotaN(1,1:3);
    ligand_anchor3_S2(1,4:10)=ligand_anchor3(1,4:10);
    
    ligand_S3(lsb1,:)=ligand_anchor2_S2;
        ligand_S3(ls,:)=ligand_anchor3_S2;
    
    
    %if ligand_S2(ls,7)~= ligand_anchor2_S22(1,7)
    %    ligand_S3(lsb1,:)=ligand_anchor2_S2;
    %    ligand_S3(ls,:)=ligand_anchor3_S2;
    %elseif ligand_S2(ls,7)== ligand_anchor2_S22(1,7) && ls < size(ligand_S2,1)
    %    ligand_S3(ls,:)=ligand_anchor2_S2;
    %    ligand_S3(ls+1,:)=ligand_anchor3_S2;
    %elseif ligand_S2(ls,7)== ligand_anchor2_S22(1,7) && ls == size(ligand_S2,1)
    %    ligand_S3(ls,:)=ligand_anchor2_S2;
    %    ligand_S3(ls-1,:)=ligand_anchor3_S2;
        
        
    %end
    elseif ligand_S2(ls,7)== ligand_anchor2_S22(1,7) && ls < size(ligand_S2,1)
        %ligand_anchor3 = ligand_S2(ls+1,:);
        continue
    elseif ligand_S2(ls,7)== ligand_anchor2_S22(1,7) && ls == size(ligand_S2,1)
        %ligand_anchor3 = ligand_S2(ls-1,:);
        continue
    end
end