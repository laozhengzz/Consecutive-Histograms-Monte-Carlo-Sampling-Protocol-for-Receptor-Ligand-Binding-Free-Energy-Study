%function forcevec = carttorquevec(protein,ligand,Paraset,centroid)
clear all
cd('C:\Users\John Zheng\Documents\MATLAB\GARF_MT\GARF_MT\GARF_potential\MT_All_Contact\potential hessian\');
home_dir = 'F:\pdb_test\';
list = dir(home_dir);
size_list = size(list);
RcutoffPL = 6;
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;

const=halva*2;
timeop=10^-12;

Nacon = 6.02*10^23;
Katommass = 1.66053886*10^-27;
for zh=3:size_list(1,1)
    tic
    if length(list(zh,1).name)<8
        continue
    elseif length(list(zh,1).name)>=8
        list_name=list(zh,1).name;
        %list_name{zh,1}=list_name;
        if strcmp(list_name(end-3:end),'.pdb')==1
            protein_id = list_name(1:end-4);
            pdbid = protein_id;
            
            protein_name = list_name;
            protein_dir = [ strcat(home_dir,protein_name)];
            [protein,protein_list_num,protein_sol] = protein_define_final_sol(protein_dir);
            
        elseif strcmp(list_name(end-4:end),'.mol2')==1
            ligand_name = list_name;
            ligand_dir = [ strcat(home_dir,ligand_name)];
            [ligand,name_str]= ligand_define_sol(ligand_dir);
        end
    end
end


run GARF_paraset
run GARF_Potential
size_Paraset=size(Paraset);

cd('C:\Users\John Zheng\Documents\MATLAB\Heatmap_Docking\');
    
    [protein_starterA,protein_starterB,repulse_num] = pocket2find_PL_Dock(protein,ligand,RcutoffPL);
    [Part_Matr_com1_new,CNN1_new,VN1_new] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d); %inter Protein-ligand interactions
    
    ZE_part_com1_new = -0.5918*Part_Matr_com1_new - VN1_new*-0.5918;
    ZN_part_com1_new = ((sqrt(CNN1_new*8+1)+1)/2-4)*log(2)+((sqrt(CNN1_new*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1_new*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1_new*8+1)+1)/2-1)*log(const);
    dG_start = (ZE_part_com1_new + ZN_part_com1_new*-0.5918)/20+3;
    cd('C:\Users\John Zheng\Documents\MATLAB\GARF_MT\GARF_MT\GARF_potential\MT_All_Contact\potential hessian\');
for loop = 1:20
    
    centroid(:,1)=sum(ligand(:,6).*ligand(:,1))./sum(ligand(:,6));
    centroid(:,2)=sum(ligand(:,6).*ligand(:,2))./sum(ligand(:,6));
    centroid(:,3)=sum(ligand(:,6).*ligand(:,3))./sum(ligand(:,6));
    
    
    clear rcl Iligand angacc_ligand ang_change force_ligand_garf
    %protein(:,1:3) = protein(:,1:3)-centroid(1,1:3);
    %ligand(:,1:3) = ligand(:,1:3)-centroid(1,1:3);
    for i = 1:size(ligand,1)
        rcl(i,1)= dist(ligand(i,:),centroid);
        Iligand(i,1) = ligand(i,6).*rcl(i,1).^2;
        %Mligand(i,1) = ligand(i,6);
    end
    
    Iligand_sum = sum(Iligand);
    Mligand_sum = sum(ligand(:,6));
    
    ligand_origin = ligand;
    
    f_ind=1;
    for i = 1:size(ligand,1)
        %for i = 6:6
        protein1 = protein;
        protein1(:,1) = protein(:,1)-ligand(i,1);
        protein1(:,2) = protein(:,2)-ligand(i,2);
        protein1(:,3) = protein(:,3)-ligand(i,3);
        
        centroid1(1,1) = centroid(:,1)-ligand(i,1);
        centroid1(1,2) = centroid(:,2)-ligand(i,3);
        centroid1(1,3) = centroid(:,3)-ligand(i,3);
        clear f_garf v_garf vf_garf torque cosangle sinangle fv_garf
        for j = 1:size(protein1,1)
            
            typea=protein1(j,4);
            typeb=ligand(i,4);
            order_ind = (typea-1)*34+typeb+1;
            rab = dist(protein1(j,:),zeros(1,3));
            rct = dist(centroid1(1,:),zeros(1,3));
            f_garf(j,1) =  -(    -Paraset(order_ind,1).*Paraset(order_ind,3).*(Paraset(order_ind,2).^Paraset(order_ind,3)) .* ( rab.^(-Paraset(order_ind,3)-1)) - Paraset(order_ind,4).*Paraset(order_ind,6).*(Paraset(order_ind,5).^Paraset(order_ind,6)) .* ( rab.^(-Paraset(order_ind,6)-1))     );
            f_garf(j,1) = (f_garf(j,1)/20).*4184;
            vf_garf(j,1) = (protein1(j,1))/rab;
            vf_garf(j,2) = (protein1(j,2))/rab;
            vf_garf(j,3) = (protein1(j,3))/rab;
            
            fv_garf(j,1) = f_garf(j,1)*vf_garf(j,1);
            fv_garf(j,2) = f_garf(j,1)*vf_garf(j,2);
            fv_garf(j,3) = f_garf(j,1)*vf_garf(j,3);
            
            ff_garf(j,1) =  -(  Paraset(order_ind,1).*Paraset(order_ind,3).*(Paraset(order_ind,3)+1).*(Paraset(order_ind,2).^Paraset(order_ind,3)) .* ( rab.^(-Paraset(order_ind,3)-2)) + Paraset(order_ind,4).*Paraset(order_ind,6).*(Paraset(order_ind,6)+1).*(Paraset(order_ind,5).^Paraset(order_ind,6)) .* ( rab.^(-Paraset(order_ind,6)-2))     );
            
            ffv_garf(j,1) = ff_garf(j,1)*vf_garf(j,1);
            ffv_garf(j,2) = ff_garf(j,1)*vf_garf(j,2);
            ffv_garf(j,3) = ff_garf(j,1)*vf_garf(j,3);
            
            vc_garf(1,1) = (-centroid1(1,1))/rct;
            vc_garf(1,2) = (-centroid1(1,2))/rct;
            vc_garf(1,3) = (-centroid1(1,3))/rct;
            
            if vf_garf(j,1) ~=0 && norm(vc_garf(1,:))~=0
                cosangle(j,1) = (vf_garf(j,1)*vc_garf(1,1) + 0*vc_garf(1,2) + 0*vc_garf(1,3)) / (abs(vf_garf(j,1))*norm(vc_garf(1,:)));
            else
                cosangle(j,1) = 1;
            end
            
            if vf_garf(j,2) ~=0 && norm(vc_garf(1,:))~=0
                cosangle(j,2) = (0*vc_garf(1,1) + vf_garf(j,2)*vc_garf(1,2) + 0*vc_garf(1,3)) / (abs(vf_garf(j,2))*norm(vc_garf(1,:)));
            else
                cosangle(j,2) = 1;
            end
            if vf_garf(j,3) ~=0 && norm(vc_garf(1,:))~=0
                cosangle(j,3) = (0*vc_garf(1,1) + 0*vc_garf(1,2) + vf_garf(j,3)*vc_garf(1,3)) / (abs(vf_garf(j,3))*norm(vc_garf(1,:)));
            else
                cosangle(j,3) = 1;
            end
            
            
            sinangle(j,1) = sqrt(1-cosangle(j,1)^2);
            sinangle(j,2) = sqrt(1-cosangle(j,2)^2);
            sinangle(j,3) = sqrt(1-cosangle(j,3)^2);
            
            
            
            %if norm(vf_garf(j,:))^2+norm(vc_garf(1,:))^2>=(dist(vf_garf(j,:),vc_garf))^2
            %sinangle = sqrt(1-cosangle^2);
            %elseif norm(vf_garf(j,:))^2+norm(vc_garf(1,:))^2<(dist(vf_garf(j,:),vc_garf))^2
            %    sinangle = sqrt(1-cosangle^2);
            %end
            
            torque(j,1) = fv_garf(j,1).*sinangle(j,1).*norm(centroid1(1,:));
            torque(j,2) = fv_garf(j,2).*sinangle(j,2).*norm(centroid1(1,:));
            torque(j,3) = fv_garf(j,3).*sinangle(j,3).*norm(centroid1(1,:));
            
            
        end
        
        torque_sum = sum(torque);
        
        angacc_ligand(i,1) = torque_sum(1,1)./(Iligand_sum*Nacon*Katommass*10^-20);
        angacc_ligand(i,2) = torque_sum(1,2)./(Iligand_sum*Nacon*Katommass*10^-20);
        angacc_ligand(i,3) = torque_sum(1,3)./(Iligand_sum*Nacon*Katommass*10^-20);
        
        ang_change(i,1) = 0.5.*angacc_ligand(i,1).*timeop^2;
        ang_change(i,2) = 0.5.*angacc_ligand(i,2).*timeop^2;
        ang_change(i,3) = 0.5.*angacc_ligand(i,3).*timeop^2;
        
        
        force_ligand_garf(i,1) = sum(fv_garf(:,1));
        force_ligand_garf(i,2) = sum(fv_garf(:,2));
        force_ligand_garf(i,3) = sum(fv_garf(:,3));
        
        
    end
    
    angacc_ligand_sum = sum(angacc_ligand);
    
    
    ang_change_sum = sum(ang_change);
    
    force_ligand_garf_sum = sum(force_ligand_garf);
    
    acc_ligand_sum = force_ligand_garf_sum./(Mligand_sum.*Nacon.*Katommass.*10^-10);
    
    dist_change_sum = 0.5.*acc_ligand_sum.*(timeop^2)./(10^-10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ang_change_sum_total(loop,:) = ang_change_sum;
    dist_change_sum_total(loop,:) = dist_change_sum;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ligand_cent_temp = ligand;
    
    ligand_cent_temp(:,1) = ligand(:,1) - centroid(1,1);
    ligand_cent_temp(:,2) = ligand(:,2) - centroid(1,2);
    ligand_cent_temp(:,3) = ligand(:,3) - centroid(1,3);
    
    ligand_rot_new_1 = rotstr(ligand_cent_temp,1,-ang_change_sum(1,1));
    ligand_rot_new_2 = rotstr(ligand_rot_new_1,2,-ang_change_sum(1,2));
    ligand_rot_new_3 = rotstr(ligand_rot_new_2,3,-ang_change_sum(1,3));
    
    ligand_rot_new = ligand_rot_new_3;
    
    ligand_rot_new(:,1) = ligand_rot_new_3(:,1) + centroid(1,1);
    ligand_rot_new(:,2) = ligand_rot_new_3(:,2) + centroid(1,2);
    ligand_rot_new(:,3) = ligand_rot_new_3(:,3) + centroid(1,3);
    
    ligand_tran_new = ligand_rot_new;
    
    ligand_tran_new(:,1) = ligand_tran_new(:,1) - dist_change_sum(1,1);
    ligand_tran_new(:,2) = ligand_tran_new(:,2) - dist_change_sum(1,2);
    ligand_tran_new(:,3) = ligand_tran_new(:,3) - dist_change_sum(1,3);
    
    cd('C:\Users\John Zheng\Documents\MATLAB\Heatmap_Docking\');
    
    [protein_starterA,protein_starterB,repulse_num] = pocket2find_PL_Dock(protein,ligand_tran_new,RcutoffPL);
    [Part_Matr_com1_new,CNN1_new,VN1_new] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d); %inter Protein-ligand interactions
    
    ZE_part_com1_new = -0.5918*Part_Matr_com1_new - VN1_new*-0.5918;
    ZN_part_com1_new = ((sqrt(CNN1_new*8+1)+1)/2-4)*log(2)+((sqrt(CNN1_new*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1_new*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1_new*8+1)+1)/2-1)*log(const);
    ligand_new_all(loop).dG = (ZE_part_com1_new + ZN_part_com1_new*-0.5918)/20+3;
    
    %if loop>1 && ligand_new_all(loop).dG>ligand_new_all(loop-1).dG
    %    dG_final = ligand_new_all(loop-1).dG;
    %    ligand_final = ligand;
    %    break
    %else
        
        ligand = ligand_tran_new;
        cd('C:\Users\John Zheng\Documents\MATLAB\GARF_MT\GARF_MT\GARF_potential\MT_All_Contact\potential hessian\');
        loop
    %end
end












