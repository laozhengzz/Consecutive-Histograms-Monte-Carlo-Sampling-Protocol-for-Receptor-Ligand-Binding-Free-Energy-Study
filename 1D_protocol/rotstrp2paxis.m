function [ligand_S3,rotaa] = rotstrp2paxis(axis,rot0,rota,rotd,ligand_S2)

%rot0 = grid_bin1(gb1,:);
%rota = ligand_anchor2_S22;
%rotd = grid_bin2(gb2,:);
rota_mark = rota;
rot0 = rot0(1,1:3);
rota = rota(1,1:3);
rotd = rotd(1,1:3);


%axis=1
warning off

ligand_S3 = zeros(size(ligand_S2,1),size(ligand_S2,2));
ligand_S3(:,4:10) = ligand_S2(:,4:10);


rota_axis=rota(1,axis);
rotd_axis=rotd(1,axis);
rot0_axis=rot0(1,axis);

rota(1,axis)=0;
rotd(1,axis)=0;
rot0(1,axis)=0;

rota_plain = rota;
rotd_plain = rotd;
rot0_plain = rot0;

rota_plain(:,axis)=[];
rotd_plain(:,axis)=[];
rot0_plain(:,axis)=[];


rota_plain_relative=rota_plain-rot0_plain;
rotd_plain_relative=rotd_plain-rot0_plain;
rot0_plain_relative=rot0_plain-rot0_plain;


da0=sqrt((rota(1,1)-rot0(1,1))^2+(rota(1,2)-rot0(1,2))^2+(rota(1,3)-rot0(1,3))^2);
dd0=sqrt((rotd(1,1)-rot0(1,1))^2+(rotd(1,2)-rot0(1,2))^2+(rotd(1,3)-rot0(1,3))^2);
dad=sqrt((rota(1,1)-rotd(1,1))^2+(rota(1,2)-rotd(1,2))^2+(rota(1,3)-rotd(1,3))^2);


if da0~=0 && dd0~=0 && dad~=0
    if rota_plain_relative(1,1)>0
        anglea0 = atan(rota_plain_relative(1,2)/rota_plain_relative(1,1));
    elseif rota_plain_relative(1,1)<0
        anglea0 = atan(rota_plain_relative(1,2)/rota_plain_relative(1,1))+pi;
    elseif rota_plain_relative(1,1)==0
        if rota_plain_relative(1,2) == 0
            anglea0 = 0;
        elseif rota_plain_relative(1,2) > 0
            anglea0 = pi/2;
        elseif rota_plain_relative(1,2) < 0
            anglea0 = -pi/2;
        end
    end
    
    if rotd_plain_relative(1,1)>0
        angled0 = atan(rotd_plain_relative(1,2)/rotd_plain_relative(1,1));
    elseif rotd_plain_relative(1,1)<0
        angled0 = atan(rotd_plain_relative(1,2)/rotd_plain_relative(1,1))+pi;
    elseif rotd_plain_relative(1,1)==0
        
        if rotd_plain_relative(1,2) == 0
            angled0 = 0;
        elseif rotd_plain_relative(1,2) > 0
            angled0 = pi/2;
        elseif rotd_plain_relative(1,2) < 0
            angled0 = -pi/2;
        end
        
    end
        
    
    delta_angle = angled0-anglea0;
    ligand_S2_coor = ligand_S2(:,1:3);
    
    ligand_S2_coor_axis=ligand_S2_coor(:,axis);
    
    ligand_S2_coor(:,axis) = 100;
    
    ligand_S2_coor_plain = ligand_S2_coor;
    
    ligand_S2_coor_plain(:,axis)=[];
    ligand_S2_coor_plain_relative(:,1)=ligand_S2_coor_plain(:,1)-rot0_plain(1,1);
    ligand_S2_coor_plain_relative(:,2)=ligand_S2_coor_plain(:,2)-rot0_plain(1,2);
    
    
    for i = 1:size(ligand_S2_coor_plain_relative,1)
        
        r_ligand_S2_coor_plain_relative = sqrt(ligand_S2_coor_plain_relative(i,1)^2+ligand_S2_coor_plain_relative(i,2)^2);
        
        if ligand_S2_coor_plain_relative(i,1)>0
            angleLS20 = atan(ligand_S2_coor_plain_relative(i,2)/ligand_S2_coor_plain_relative(i,1));
        elseif ligand_S2_coor_plain_relative(i,1)<0
            angleLS20 = atan(ligand_S2_coor_plain_relative(i,2)/ligand_S2_coor_plain_relative(i,1))+pi;
        elseif ligand_S2_coor_plain_relative(i,1)==0
            
            if ligand_S2_coor_plain_relative(i,2) == 0
                angleLS20 = 0;
            elseif ligand_S2_coor_plain_relative(i,2) > 0
                angleLS20 = pi/2;
            elseif ligand_S2_coor_plain_relative(i,2) < 0
                angleLS20 = -pi/2;
            end
            
        end
        
        angleLS20_new = angleLS20 + delta_angle;
        
        ligand_S2_coor_plain_relative_new(i,1) = r_ligand_S2_coor_plain_relative*cos(angleLS20_new);
        ligand_S2_coor_plain_relative_new(i,2) = r_ligand_S2_coor_plain_relative*sin(angleLS20_new);
        mk=1;
        for j = 1:size(ligand_S2_coor,2)
            if ligand_S2_coor(i,j)~=100 && mk == 1
                ligand_S2_coor(i,j) = ligand_S2_coor_plain_relative_new(i,1)+rot0_plain(1,1);
                mk=2;
            elseif ligand_S2_coor(i,j)~=100 && mk == 2
                ligand_S2_coor(i,j) = ligand_S2_coor_plain_relative_new(i,2)+rot0_plain(1,2);
                mk=3;
            elseif ligand_S2_coor(i,j)==100
                ligand_S2_coor(i,j) = ligand_S2_coor_axis(i,1);
                
            end
        end
        
    end
    
    
    ligand_S3(:,1:3) = ligand_S2_coor;
    
elseif da0==0 || dd0==0 || dad==0
    ligand_S3 = ligand_S2;
end


rotaa = ligand_S3(ligand_S3(:,7)==rota_mark(:,7),:);

