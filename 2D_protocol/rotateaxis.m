function [rotaN,ligand_anchor3_rotaN] = rotateaxis(axis,rot0,rota,rotd,ligand_anchor3_S1)
%axis=1
warning off
rotaN=[];
ligand_anchor3_rotaN=[];
%ligand_anchor1_S2 = ligand_anchor1_S1;
%rota=rotaN;
%rotd=rotd;
%rot0=rot0;
%ligand_anchor3_S1=ligand_anchor3_rotaN;

rota_axis=rota(1,axis);
rotd_axis=rotd(1,axis);
rot0_axis=rot0(1,axis);
ligand_anchor3_S1_axis=ligand_anchor3_S1(1,axis);

rota(1,axis)=0;
rotd(1,axis)=0;
rot0(1,axis)=0;
ligand_anchor3_S1(1,axis)=0;

da0=sqrt((rota(1,1)-rot0(1,1))^2+(rota(1,2)-rot0(1,2))^2+(rota(1,3)-rot0(1,3))^2);
dd0=sqrt((rotd(1,1)-rot0(1,1))^2+(rotd(1,2)-rot0(1,2))^2+(rotd(1,3)-rot0(1,3))^2);
dad=sqrt((rota(1,1)-rotd(1,1))^2+(rota(1,2)-rotd(1,2))^2+(rota(1,3)-rotd(1,3))^2);


if da0~=0 && dd0~=0 && dad~=0


% distance anchor3_S1 to rot0
d30=sqrt((ligand_anchor3_S1(1,1)-rot0(1,1))^2+(ligand_anchor3_S1(1,2)-rot0(1,2))^2+(ligand_anchor3_S1(1,3)-rot0(1,3))^2);
if axis==1
    rotr=[rot0(1,1),rot0(1,2)+da0,rot0(1,3)];
elseif axis==2
    rotr=[rot0(1,1)+da0,rot0(1,2),rot0(1,3)];
elseif axis==3
    rotr=[rot0(1,1)+da0,rot0(1,2),rot0(1,3)];
end

dr0=da0;
dar=sqrt((rota(1,1)-rotr(1,1))^2+(rota(1,2)-rotr(1,2))^2+(rota(1,3)-rotr(1,3))^2);
ddr=sqrt((rotd(1,1)-rotr(1,1))^2+(rotd(1,2)-rotr(1,2))^2+(rotd(1,3)-rotr(1,3))^2);
d3r=sqrt((ligand_anchor3_S1(1,1)-rotr(1,1))^2+(ligand_anchor3_S1(1,2)-rotr(1,2))^2+(ligand_anchor3_S1(1,3)-rotr(1,3))^2);

if axis==1
    
    if rota(1,3)-rot0(1,3)>=0
        anglear=acos((da0^2+dr0^2-dar^2)/(2*da0*dr0));
    elseif rota(1,3)-rot0(1,3)<0
        anglear=2*pi-acos((da0^2+dr0^2-dar^2)/(2*da0*dr0));
    end
    
    if rotd(1,3)-rot0(1,3)>=0
        angledr=acos((dd0^2+dr0^2-ddr^2)/(2*dd0*dr0));
    elseif rotd(1,3)-rot0(1,3)<0
        angledr=2*pi-acos((dd0^2+dr0^2-ddr^2)/(2*dd0*dr0));
    end
    
    if ligand_anchor3_S1(1,3)-rot0(1,3)>=0
        angle3r=acos((d30^2+dr0^2-d3r^2)/(2*d30*dr0));
    elseif ligand_anchor3_S1(1,3)-rot0(1,3)<0
        angle3r=2*pi-acos((d30^2+dr0^2-d3r^2)/(2*d30*dr0));
    end
    
elseif axis==2
    
    if rota(1,3)-rot0(1,3)>=0
        anglear=acos((da0^2+dr0^2-dar^2)/(2*da0*dr0));
    elseif rota(1,3)-rot0(1,3)<0
        anglear=2*pi-acos((da0^2+dr0^2-dar^2)/(2*da0*dr0));
    end
    
    if rotd(1,3)-rot0(1,3)>=0
        angledr=acos((dd0^2+dr0^2-ddr^2)/(2*dd0*dr0));
    elseif rotd(1,3)-rot0(1,3)<0
        angledr=2*pi-acos((dd0^2+dr0^2-ddr^2)/(2*dd0*dr0));
    end
    
    if ligand_anchor3_S1(1,3)-rot0(1,3)>=0
        angle3r=acos((d30^2+dr0^2-d3r^2)/(2*d30*dr0));
    elseif ligand_anchor3_S1(1,3)-rot0(1,3)<0
        angle3r=2*pi-acos((d30^2+dr0^2-d3r^2)/(2*d30*dr0));
    end
    
elseif axis==3
    
    if rota(1,2)-rot0(1,2)>=0
        anglear=acos((da0^2+dr0^2-dar^2)/(2*da0*dr0));
    elseif rota(1,2)-rot0(1,2)<0
        anglear=2*pi-acos((da0^2+dr0^2-dar^2)/(2*da0*dr0));
    end
    
    if rotd(1,2)-rot0(1,2)>=0
        angledr=acos((dd0^2+dr0^2-ddr^2)/(2*dd0*dr0));
    elseif rotd(1,2)-rot0(1,2)<0
        angledr=2*pi-acos((dd0^2+dr0^2-ddr^2)/(2*dd0*dr0));
    end
    
    if ligand_anchor3_S1(1,2)-rot0(1,2)>=0
        angle3r=acos((d30^2+dr0^2-d3r^2)/(2*d30*dr0));
    elseif ligand_anchor3_S1(1,2)-rot0(1,2)<0
        angle3r=2*pi-acos((d30^2+dr0^2-d3r^2)/(2*d30*dr0));
    end
    
end
anglead=angledr-anglear;
anglead = real(anglead);
angleanf = anglear+anglead;
angle3nf = angle3r+anglead;

angleanf = real(angleanf);
angle3nf = real(angle3nf);

if axis==1
    if (da0*cos(angleanf))/(rotd(1,2)-rot0(1,2))>=0
        rotaN(1,2)=da0*cos(angleanf)+rot0(1,2);
        ligand_anchor3_rotaN(1,2)=d30*cos(angle3nf)+rot0(1,2);
    elseif (da0*cos(angleanf))/(rotd(1,2)-rot0(1,2))<0
        rotaN(1,2)=-da0*cos(angleanf)+rot0(1,2);
        ligand_anchor3_rotaN(1,2)=-d30*cos(angle3nf)+rot0(1,2);
    elseif abs(rotd(1,2)-rot0(1,2))<10^-7
        rotaN(1,2)=da0*cos(angleanf)+rot0(1,2);
        ligand_anchor3_rotaN(1,2)=d30*cos(angle3nf)+rot0(1,2);
    end
    
    if (da0*sin(angleanf))/(rotd(1,3)-rot0(1,3))>=0
        rotaN(1,3)=da0*sin(angleanf)+rot0(1,3);
        ligand_anchor3_rotaN(1,3)=d30*sin(angle3nf)+rot0(1,3);
    elseif (da0*sin(angleanf))/(rotd(1,3)-rot0(1,3))<0
        rotaN(1,3)=-da0*sin(angleanf)+rot0(1,3);
        ligand_anchor3_rotaN(1,3)=-d30*sin(angle3nf)+rot0(1,3);
    elseif abs(rotd(1,3)-rot0(1,3))<10^-7    
        rotaN(1,3)=da0*sin(angleanf)+rot0(1,3);
        ligand_anchor3_rotaN(1,3)=d30*sin(angle3nf)+rot0(1,3);
        
    end
elseif axis==2
    if (da0*cos(angleanf))/(rotd(1,1)-rot0(1,1))>=0
        rotaN(1,1)=da0*cos(angleanf)+rot0(1,1);
        ligand_anchor3_rotaN(1,1)=d30*cos(angle3nf)+rot0(1,1);
    elseif (da0*cos(angleanf))/(rotd(1,1)-rot0(1,1))<0
        rotaN(1,1)=-da0*cos(angleanf)+rot0(1,1);
        ligand_anchor3_rotaN(1,1)=-d30*cos(angle3nf)+rot0(1,1);
    elseif abs(rotd(1,1)-rot0(1,1))<10^-7
        rotaN(1,1)=da0*cos(angleanf)+rot0(1,1);
        ligand_anchor3_rotaN(1,1)=d30*cos(angle3nf)+rot0(1,1);
        
    end
    
    if (da0*sin(angleanf))/(rotd(1,3)-rot0(1,3))>=0
        rotaN(1,3)=da0*sin(angleanf)+rot0(1,3);
        ligand_anchor3_rotaN(1,3)=d30*sin(angle3nf)+rot0(1,3);
    elseif (da0*sin(angleanf))/(rotd(1,3)-rot0(1,3))<0
        rotaN(1,3)=-da0*sin(angleanf)+rot0(1,3);
        ligand_anchor3_rotaN(1,3)=-d30*sin(angle3nf)+rot0(1,3);
    elseif abs(rotd(1,3)-rot0(1,3))<10^-7
        rotaN(1,3)=da0*sin(angleanf)+rot0(1,3);
        ligand_anchor3_rotaN(1,3)=d30*sin(angle3nf)+rot0(1,3);
        
    end
    
elseif axis==3
    if (da0*cos(angleanf))/(rotd(1,1)-rot0(1,1))>=0
        rotaN(1,1)=da0*cos(angleanf)+rot0(1,1);
        ligand_anchor3_rotaN(1,1)=d30*cos(angle3nf)+rot0(1,1);
    elseif (da0*cos(angleanf))/(rotd(1,1)-rot0(1,1))<0
        rotaN(1,1)=-da0*cos(angleanf)+rot0(1,1);
        ligand_anchor3_rotaN(1,1)=-d30*cos(angle3nf)+rot0(1,1);
        
    elseif abs(rotd(1,1)-rot0(1,1))<10^-7
        rotaN(1,1)=da0*cos(angleanf)+rot0(1,1);
        ligand_anchor3_rotaN(1,1)=d30*cos(angle3nf)+rot0(1,1);
        
    end
    
    if (da0*sin(angleanf))/(rotd(1,2)-rot0(1,2))>=0
        rotaN(1,2)=da0*sin(angleanf)+rot0(1,2);
        ligand_anchor3_rotaN(1,2)=d30*sin(angle3nf)+rot0(1,2);
    elseif (da0*sin(angleanf))/(rotd(1,2)-rot0(1,2))<0
        rotaN(1,2)=-da0*sin(angleanf)+rot0(1,2);
        ligand_anchor3_rotaN(1,2)=-d30*sin(angle3nf)+rot0(1,2);
        
    elseif abs(rotd(1,2)-rot0(1,2))<10^-7
        rotaN(1,2)=da0*sin(angleanf)+rot0(1,2);
        ligand_anchor3_rotaN(1,2)=d30*sin(angle3nf)+rot0(1,2);
        
    end
end


rotaN(1,axis)=rota_axis;
ligand_anchor3_rotaN(1,axis)=ligand_anchor3_S1_axis;


elseif da0==0 || dd0==0 || dad==0
    rotaN = rota;
    rotaN(1,axis)=rota_axis;
    ligand_anchor3_rotaN = ligand_anchor3_S1;
    ligand_anchor3_rotaN(1,axis)=ligand_anchor3_S1_axis;
end
    


