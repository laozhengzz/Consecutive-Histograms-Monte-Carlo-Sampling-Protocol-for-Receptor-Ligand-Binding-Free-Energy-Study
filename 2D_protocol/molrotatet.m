function [ligand_S22,ligand_anchor2_S22] = molrotatet(ligand_S2,ligand_anchor1_S2,ligand_anchor2_S2,anglet)       
                                           
%ligand_S2 = ligand_S22;
ligand_S2(:,1) = ligand_S2(:,1) + (0-ligand_anchor1_S2(1,1));
ligand_S2(:,2) = ligand_S2(:,2) + (0-ligand_anchor1_S2(1,2));
ligand_S2(:,3) = ligand_S2(:,3) + (0-ligand_anchor1_S2(1,3));

for i = 1:size(ligand_S2,1)
    
    
    rcoor = sqrt(ligand_S2(i,1)^2+ligand_S2(i,2)^2+ligand_S2(i,3)^2);
    theta = acos(   ligand_S2(i,3)/sqrt(ligand_S2(i,1)^2+ligand_S2(i,2)^2+ligand_S2(i,3)^2) );
    if ligand_S2(i,1)>0
        phi = atan(ligand_S2(i,2)/ligand_S2(i,1));
    elseif ligand_S2(i,1)<0
        phi = atan(ligand_S2(i,2)/ligand_S2(i,1))+pi;
    elseif ligand_S2(i,1)==0
        
        if ligand_S2(i,2)>0
            phi = pi/2;
        elseif ligand_S2(i,2)<0
            phi = pi/2+pi;
        elseif ligand_S2(i,2) == 0
            phi = 0;
        end
        
    end
    
    phi = phi+anglet;
    
    if phi~=0
        
        ligand_S22(i,1) = rcoor*sin(theta)*cos(phi);
        ligand_S22(i,2) = rcoor*sin(theta)*sin(phi);
        ligand_S22(i,3) = rcoor*cos(theta);
        
        
        ligand_S22(i,1) = ligand_S22(i,1) - (0-ligand_anchor1_S2(1,1));
        ligand_S22(i,2) = ligand_S22(i,2) - (0-ligand_anchor1_S2(1,2));
        ligand_S22(i,3) = ligand_S22(i,3) - (0-ligand_anchor1_S2(1,3));
        
        ligand_S22(i,4:10) = ligand_S2(i,4:10);
        
    elseif phi==0
        ligand_S22 = ligand_S2;
        ligand_S22(i,1) = ligand_S22(i,1) - (0-ligand_anchor1_S2(1,1));
        ligand_S22(i,2) = ligand_S22(i,2) - (0-ligand_anchor1_S2(1,2));
        ligand_S22(i,3) = ligand_S22(i,3) - (0-ligand_anchor1_S2(1,3));
        
        ligand_S22(i,4:10) = ligand_S2(i,4:10);
    end
    
end
    
ligand_anchor2_S22=ligand_S22(ligand_S22(:,7)==ligand_anchor2_S2(1,7),:);
    
