function ligand_output = cartvecNRG_waterF(protein,ligand,loop_num)
%cd('/home/dell/Octave_work/GARF_MT/GARF_MT/GARF_potential/MT_All_Contact/potential_hessian/');

RcutoffPL = 6;
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;

const=halva*2;
timeop=10^-13;
NR_step = 10^-0;

Nacon = 6.02*10^23;
Katommass = 1.66053886*10^-27;

ligand_ori = ligand;
run GARF_paraset
run GARF_Potential
size_Paraset=size(Paraset);

%cd('/home/dell/Octave_work/Heatmap_Docking/');
    
    [protein_starterA,protein_starterB,repulse_num] = pocket2find_PL_Dock_water(protein,ligand,RcutoffPL);
    [Part_Matr_com1_new,CNN1_new,VN1_new] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d); %inter Protein-ligand interactions
    
    ZE_part_com1_new = -0.5918*Part_Matr_com1_new - VN1_new*-0.5918;
    ZN_part_com1_new = ((sqrt(CNN1_new*8+1)+1)/2-4)*log(2)+((sqrt(CNN1_new*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1_new*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1_new*8+1)+1)/2-1)*log(const);
    dG_start = (ZE_part_com1_new + ZN_part_com1_new*-0.5918)/20+3;
%    cd('/home/dell/Octave_work/GARF_MT/GARF_MT/GARF_potential/MT_All_Contact/potential_hessian/');
for loop = 1:loop_num
    
    centroid(:,1)=sum(ligand(:,6).*ligand(:,1))./sum(ligand(:,6));
    centroid(:,2)=sum(ligand(:,6).*ligand(:,2))./sum(ligand(:,6));
    centroid(:,3)=sum(ligand(:,6).*ligand(:,3))./sum(ligand(:,6));
    
    
    clear rcl Iligand angacc_ligand ang_change force_ligand_garf fv_garf_LAtom torque_LAtom ffv_garf_LAtom ftorque_LAtom
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
    ffv_garf_L_final = zeros(3,3);
    affv_garf_L_final = zeros(3,3);
    ftorque_L_final = zeros(3,3);
    aftorque_L_final = zeros(3,3);
    
    f_ind=1;
    for i = 1:size(ligand,1)
    %for i = 9:9
        protein1 = protein;
        protein1(:,1) = protein(:,1)-ligand(i,1);
        protein1(:,2) = protein(:,2)-ligand(i,2);
        protein1(:,3) = protein(:,3)-ligand(i,3);
        
        centroid1(1,1) = centroid(:,1)-ligand(i,1);
        centroid1(1,2) = centroid(:,2)-ligand(i,3);
        centroid1(1,3) = centroid(:,3)-ligand(i,3);
        clear f_garf v_garf vf_garf vff_garf torque cosangle sinangle fv_garf ffv_garf_L PL_Potential
        
        ffv_garf_L = zeros(3,3);
        affv_garf_L = zeros(3,3);
        
        ftorque_L = zeros(3,3);
        faang_L = zeros(3,3);
        
        for j = 1:size(protein1,1)
            
            typea=protein1(j,4);
            typeb=ligand(i,4);
            order_ind = (typea-1)*34+typeb+1;
            rab = dist(protein1(j,:),zeros(1,3));
            rct = dist(centroid1(1,:),zeros(1,3));
            vf_garf(j,1) = (protein1(j,1))/rab;
            vf_garf(j,2) = (protein1(j,2))/rab;
            vf_garf(j,3) = (protein1(j,3))/rab;
            
            vff_garf(1,1) = 1/rab - (protein1(j,1)^2)/(rab^3);
            vff_garf(1,2) = - protein1(j,1).*protein1(j,2)/(rab^3);
            vff_garf(1,3) = - protein1(j,1).*protein1(j,3)/(rab^3);
            
            vff_garf(2,1) = - protein1(j,2).*protein1(j,1)/(rab^3);
            vff_garf(2,2) = 1/rab - (protein1(j,2)^2)/(rab^3);
            vff_garf(2,3) = - protein1(j,2).*protein1(j,3)/(rab^3);
            
            vff_garf(3,1) = - protein1(j,3).*protein1(j,1)/(rab^3);
            vff_garf(3,2) = - protein1(j,3).*protein1(j,2)/(rab^3);
            vff_garf(3,3) = 1/rab - (protein1(j,3)^2)/(rab^3);
            
            PL_Potential(j,1) = (Paraset(order_ind,1).*(Paraset(order_ind,2)./rab).^Paraset(order_ind,3)+Paraset(order_ind,4).*(Paraset(order_ind,5)./rab).^Paraset(order_ind,6));
            
            for fi = 1:3
                
                
                
                fv_garf(j,fi) =  -(    -Paraset(order_ind,1).*Paraset(order_ind,3).*(Paraset(order_ind,2).^Paraset(order_ind,3)) .* ( rab.^(-Paraset(order_ind,3)-1)) - Paraset(order_ind,4).*Paraset(order_ind,6).*(Paraset(order_ind,5).^Paraset(order_ind,6)) .* ( rab.^(-Paraset(order_ind,6)-1))     ).*vf_garf(j,fi);
                fv_garf(j,fi) = (fv_garf(j,fi)/20).*4184;
                
                afv_garf(j,fi) = fv_garf(j,fi)./(ligand(i,6).*Nacon.*Katommass);
                
            end
                
            for fi = 1:3
                for fii = 1:3
                
                    ffv_garf(fi,fii) =  -(  (Paraset(order_ind,1).*Paraset(order_ind,3).*(Paraset(order_ind,3)+1).*(Paraset(order_ind,2).^Paraset(order_ind,3)) .* ( rab.^(-Paraset(order_ind,3)-2)) + Paraset(order_ind,4).*Paraset(order_ind,6).*(Paraset(order_ind,6)+1).*(Paraset(order_ind,5).^Paraset(order_ind,6)) .* ( rab.^(-Paraset(order_ind,6)-2))).*vf_garf(j,fii).*vf_garf(j,fi)      ...
                        + (    -Paraset(order_ind,1).*Paraset(order_ind,3).*(Paraset(order_ind,2).^Paraset(order_ind,3)) .* ( rab.^(-Paraset(order_ind,3)-1)) - Paraset(order_ind,4).*Paraset(order_ind,6).*(Paraset(order_ind,5).^Paraset(order_ind,6)) .* ( rab.^(-Paraset(order_ind,6)-1))     ).*vff_garf(fi,fii)  );
                    ffv_garf(fi,fii) = (ffv_garf(fi,fii)/20).*4184;
                    
                    affv_garf(fi,fii) = ffv_garf(fi,fii)./(ligand(i,6).*Nacon.*Katommass);
                    
                end
            end
                
            ffv_garf_L = ffv_garf_L + ffv_garf;
            affv_garf_L = affv_garf_L + affv_garf;
                %fv_garf(j,fi) = real(fv_garf(j,fi));
                
                
                
                %ffv_garf(j,fi) = real(ffv_garf(j,fi));
                
            %vc_garf(1,1) = (-centroid1(1,1))/rct;
            %vc_garf(1,2) = (-centroid1(1,2))/rct;
            %vc_garf(1,3) = (-centroid1(1,3))/rct;
            
            torque(j,1) = fv_garf(j,2).*-(centroid1(1,3)) + fv_garf(j,3).*-(centroid1(1,2));
            torque(j,2) = fv_garf(j,1).*-(centroid1(1,3)) + fv_garf(j,3).*-(centroid1(1,1));
            torque(j,3) = fv_garf(j,1).*-(centroid1(1,2)) + fv_garf(j,2).*-(centroid1(1,1));
            
            aang(j,1) = (fv_garf(j,2).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (fv_garf(j,3).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2);
            aang(j,2) = (fv_garf(j,1).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (fv_garf(j,3).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            aang(j,3) = (fv_garf(j,1).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2) + (fv_garf(j,2).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            
            
            ftorque(1,1) = ffv_garf(2,1).*-(centroid1(1,3)) + ffv_garf(3,1).*-(centroid1(1,2));
            ftorque(1,2) = ffv_garf(2,2).*-(centroid1(1,3)) + ffv_garf(3,2).*-(centroid1(1,2));
            ftorque(1,3) = ffv_garf(2,3).*-(centroid1(1,3)) + ffv_garf(3,3).*-(centroid1(1,2));
            
            ftorque(2,1) = ffv_garf(1,1).*-(centroid1(1,3)) + ffv_garf(3,1).*-(centroid1(1,1));
            ftorque(2,2) = ffv_garf(1,2).*-(centroid1(1,3)) + ffv_garf(3,2).*-(centroid1(1,1));
            ftorque(2,3) = ffv_garf(1,3).*-(centroid1(1,3)) + ffv_garf(3,3).*-(centroid1(1,1));
            
            ftorque(3,1) = ffv_garf(1,1).*-(centroid1(1,2)) + ffv_garf(2,1).*-(centroid1(1,1));
            ftorque(3,2) = ffv_garf(1,2).*-(centroid1(1,2)) + ffv_garf(2,2).*-(centroid1(1,1));
            ftorque(3,3) = ffv_garf(1,3).*-(centroid1(1,2)) + ffv_garf(2,3).*-(centroid1(1,1));
            
            
            ftorque_L = ftorque_L + ftorque;
            
            
            
            faang(1,1) = (ffv_garf(2,1).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (ffv_garf(3,1).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2);
            faang(1,2) = (ffv_garf(2,2).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (ffv_garf(3,2).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2);
            faang(1,3) = (ffv_garf(2,3).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (ffv_garf(3,3).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2);
            
            faang(2,1) = (ffv_garf(1,1).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (ffv_garf(3,1).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            faang(2,2) = (ffv_garf(1,2).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (ffv_garf(3,2).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            faang(2,3) = (ffv_garf(1,3).*-(centroid1(1,3)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,3))^2) + (ffv_garf(3,3).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            
            faang(3,1) = (ffv_garf(1,1).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2) + (ffv_garf(2,1).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            faang(3,2) = (ffv_garf(1,2).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2) + (ffv_garf(2,2).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            faang(3,3) = (ffv_garf(1,3).*-(centroid1(1,2)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,2))^2) + (ffv_garf(2,3).*-(centroid1(1,1)))./(ligand(i,6).*Nacon.*Katommass.*(centroid1(1,1))^2);
            
            faang_L = faang_L + faang;
            
        end
        
        fv_garf_LAtom(i,:) = sum(fv_garf);
        torque_LAtom(i,:) = sum(torque);
        
        afv_garf_LAtom(i,:) = sum(afv_garf);
        aang_LAtom(i,:) = sum(aang);
        
        
        aang_LAtom2(i,:) = aang_LAtom(i,:)*(faang_L^-1);
        afv_garf_LAtom2(i,:) = afv_garf_LAtom(i,:)*(affv_garf_L^-1);
        
        
        ang_change(i,1) = 0.5.*aang_LAtom(i,1).*timeop^2;
        ang_change(i,2) = 0.5.*aang_LAtom(i,2).*timeop^2;
        ang_change(i,3) = 0.5.*aang_LAtom(i,3).*timeop^2;
        
        dist_change(i,1) = 0.5*afv_garf_LAtom(i,1).*timeop^2;
        dist_change(i,2) = 0.5*afv_garf_LAtom(i,2).*timeop^2;
        dist_change(i,3) = 0.5*afv_garf_LAtom(i,3).*timeop^2;
        
        %[eigVecff,eigValff]=eig(ffv_garf_L);
        %[eigVecft,eigValft]=eig(ftorque_L);
        
        %dist_change(i,:) = ((ffv_garf_L^-1)*(fv_garf_LAtom(i,:)'))';
        %dist_change1(i,:) = fv_garf_LAtom(i,:);
        
        %ang_change(i,:) = ((ftorque_L^-1)*(torque_LAtom(i,:)'))';
        %ang_change1(i,:) = torque_LAtom(i,:);
        
        PL_Potential_L(i,1)=sum(PL_Potential);
        
        %ffv_garf_L_final(i).ffv_garf_L = ffv_garf_L;
        
        %ftorque_L_final = ftorque_L_final + ftorque_L;
    end
    
    
    fv_garf_final = sum(fv_garf_LAtom);
    torque_final = sum(torque_LAtom);
    
    
    afv_garf_final = sum(afv_garf_LAtom);
    aang_final = sum(aang_LAtom);
    
    dist_change_ff = sum(dist_change);
    ang_change_ff = sum(ang_change);
    
    
    
    for ffi = 1:40
        fpara = 10^ffi;
        dist_change_ftemp=dist_change_ff.*fpara;
        ang_change_ftemp=ang_change_ff.*fpara;
        if (abs(dist_change_ftemp(1,1))>0.1 && abs(dist_change_ftemp(1,1))<1) || (abs(dist_change_ftemp(1,2))>0.1 && abs(dist_change_ftemp(1,2))<1) || (abs(dist_change_ftemp(1,3))>0.1 && abs(dist_change_ftemp(1,3))<1) || (abs(ang_change_ftemp(1,1))>0.1 && abs(ang_change_ftemp(1,1))<1) || (abs(ang_change_ftemp(1,2))>0.1 && abs(ang_change_ftemp(1,2))<1) || (abs(ang_change_ftemp(1,3))>0.1 && abs(ang_change_ftemp(1,3))<1)
            dist_change_ff = dist_change_ftemp;
            ang_change_ff = ang_change_ftemp;
            break
        end
    end
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ri_end = 5;
    ri_ind = 1;
    %ri = 1;
    clear ligand_new_ff ligand_new_ff_NR_dG
    %while ri <= ri_end || (ri_ind == 0 && ri <= 10)
    for ri = -2:ri_end
        rpara = (0.01)*10^(ri-1);
        dist_change_final = dist_change_ff*rpara;
        ang_change_final = ang_change_ff*rpara;
        
        
        NR_ang_change_sum = ang_change_final;
        
        NR_dist_change_sum = dist_change_final;
        
        
        %NR_ang_change_sum = NR_step.*fv_garf_sum_norm./para_pwer.*ang_change_sum./fang_change_sum;
        
        %NR_dist_change_sum = NR_step.*fv_garf_sum_norm./para_pwer.*force_ligand_garf_sum./fforce_ligand_garf_sum;
        
        %NR_ang_change_sum = NR_step.*fv_garf_sum_norm.*ang_change_sum./abs(fang_change_sum);
        
        %NR_dist_change_sum = NR_step.*fv_garf_sum_norm.*dist_change_sum./abs(fdist_change_sum);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ligand_cent_temp = ligand;
        
        ligand_cent_temp(:,1) = ligand(:,1) - centroid(1,1);
        ligand_cent_temp(:,2) = ligand(:,2) - centroid(1,2);
        ligand_cent_temp(:,3) = ligand(:,3) - centroid(1,3);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        NR_ligand_rot_new_1 = rotstr(ligand_cent_temp,1,-NR_ang_change_sum(1,1));
        NR_ligand_rot_new_2 = rotstr(NR_ligand_rot_new_1,2,-NR_ang_change_sum(1,2));
        NR_ligand_rot_new_3 = rotstr(NR_ligand_rot_new_2,3,-NR_ang_change_sum(1,3));
        
        NR_ligand_rot_new = NR_ligand_rot_new_3;
        
        NR_ligand_rot_new(:,1) = NR_ligand_rot_new_3(:,1) + centroid(1,1);
        NR_ligand_rot_new(:,2) = NR_ligand_rot_new_3(:,2) + centroid(1,2);
        NR_ligand_rot_new(:,3) = NR_ligand_rot_new_3(:,3) + centroid(1,3);
        
        NR_ligand_tran_new = NR_ligand_rot_new;
        
        NR_ligand_tran_new(:,1) = NR_ligand_tran_new(:,1) - NR_dist_change_sum(1,1);
        NR_ligand_tran_new(:,2) = NR_ligand_tran_new(:,2) - NR_dist_change_sum(1,2);
        NR_ligand_tran_new(:,3) = NR_ligand_tran_new(:,3) - NR_dist_change_sum(1,3);
        
        %rmsd_NR(loop,1) = RMSD_cal(NR_ligand_tran_new,ligand_ori);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%        cd('/home/dell/Octave_work/Heatmap_Docking/');
        
        [protein_starterA,protein_starterB,repulse_num] = pocket2find_PL_Dock(protein,NR_ligand_tran_new,RcutoffPL);
        
        if ~isempty(protein_starterA) && ~isempty(protein_starterB)
            
            [NR_Part_Matr_com1_new,NR_CNN1_new,NR_VN1_new] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d); %inter Protein-ligand interactions
            
            NR_ZE_part_com1_new = -0.5918*NR_Part_Matr_com1_new - NR_VN1_new*-0.5918;
            NR_ZN_part_com1_new = ((sqrt(NR_CNN1_new*8+1)+1)/2-4)*log(2)+((sqrt(NR_CNN1_new*8+1)+1)/2-3)*log(2*pi)+((sqrt(NR_CNN1_new*8+1)+1)/2-2)*log(4*pi)+((sqrt(NR_CNN1_new*8+1)+1)/2-1)*log(const);
            NR_dG = (NR_ZE_part_com1_new + NR_ZN_part_com1_new*-0.5918)/20+3;
            
            ligand_new_ff_NR_dG(ri_ind,1) = NR_dG;
            ligand_new_ff(ri_ind).ligand = NR_ligand_tran_new;
            ri_ind = ri_ind + 1;
            %ri = ri+1;
        elseif isempty(protein_starterA) || isempty(protein_starterB)
            if ri == 1
                ligand_new_ff_NR_dG(ri_ind,1) = 100;
                ligand_new_ff(ri_ind).ligand = [];
                ri_ind = ri_ind + 1;
                %ri = ri+1;
                continue
            elseif ri>1
                ligand_new_ff_NR_dG(ri_ind,1) = 100;
                ligand_new_ff(ri_ind).ligand = [];
                ri_ind = ri_ind + 1;
                %ri = ri+1;
                continue
            end
        end
        
        
    end
    
    [ff1,ff2] = find(ligand_new_ff_NR_dG == min(ligand_new_ff_NR_dG));
    ligand_new_all_NR_dG(loop,1) = ligand_new_ff_NR_dG(ff1(1,1),1);
    ligand_new_all(loop).ligand = ligand_new_ff(ff1(1,1)).ligand;
    [repulse_num] = repulsenum_water(protein,ligand_new_ff(ff1(1,1)).ligand)
    ligand_new_all(loop).repulse_num = repulse_num;
    
    ligand = ligand_new_ff(ff1(1,1)).ligand;
        %ligand = ligand_tran_new;
        %ligand = NR_ligand_tran_new;
%        cd('/home/dell/Octave_work/GARF_MT/GARF_MT/GARF_potential/MT_All_Contact/potential_hessian/');
        loop;
    %end
    
end


[ml1,ml2] = find(ligand_new_all_NR_dG == min(ligand_new_all_NR_dG));

ligand_output.dG = ligand_new_all_NR_dG(ml1,1);
ligand_output.coordinate = ligand_new_all(ml1).ligand;
ligand_output.repulse_num = ligand_new_all(ml1).repulse_num;











