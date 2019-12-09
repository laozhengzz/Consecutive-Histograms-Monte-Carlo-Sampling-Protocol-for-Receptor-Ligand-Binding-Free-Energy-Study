clear all

home_dir = strcat('F:\pdbbind_v2018\pdbbind_v2018\pdbbind_v2018_refined\refined-set\1bty\');
list1=dir(home_dir);
size_list1=size(list1);
%Aout=zeros(200,6);
%Bout=zeros(200,28);
dist_cent_search = 4;
RcutoffPL=6;
R_srch=0.5;
MTsample = R_srch;
R_L_tor_srch=0.1;
rd = 0.005;
hnstate = R_srch/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2;


Rcutoff=5;
tn=0;
clear test_name
%for zh=15:15
code_dir = pwd;
%jpz=1;
%for zh=floor(size_list(1,1)*5/6)+1:size_list(1,1)
run GARF_Potential

for zh=1:size_list1(1,1)
    tic
    if length(list1(zh,1).name)<8
        continue
    elseif length(list1(zh,1).name)>8
        list1_name=list1(zh,1).name;
        %list_name{zh,1}=list_name;
        if strcmp(list1_name(end-3:end),'.pdb')==1 %&& strcmp(list1_name(1:4),list_name)==1
            protein_id = list1_name(1:4);
            pdbid = protein_id;
            
            protein_name = list1_name;
            protein_dir = [ strcat(home_dir,list1_name)];
            [protein,protein_Namelist] = protein_define_nowater(protein_dir);
        end
    end
end

%39
%for zh=1:size_list1(1,1)
for zh=1:size_list1(1,1)
    %for zh=5:25
    tic
    if length(list1(zh,1).name)<8
        continue
    elseif length(list1(zh,1).name)>8
        list1_name=list1(zh,1).name;
        
        if strcmp(list1_name(end-4:end),'.mol2')==1
            ligand_name = [ strcat(list1_name)];
            ligand_dir = [ strcat(home_dir,list1_name)];
            
            
            
            %[pocket,pocket_Namelist] = protein_define_final(pocket_dir);
            [ligand,name_str]= ligand_define_sol_spec(ligand_dir);
            ligand_origin = ligand;
            %ligand(ligand(:,4)==6,4)=5;
            
            %[ZE_P_san] = KMTISM_function_code(protein);
            %[ZE_L_san] = KMTISM_function_code(ligand);
            
            %boundary_id = [ strcat(ligand_name(1:end-5),'-blob')];
            %boundary_name = [ strcat(boundary_id,'.txt')];
            %boundary_dir = [ strcat(home_dir,pdbid,'\',boundary_name)];
            
            %[region]= boundary_define_simp(boundary_dir);
            
            %if ligpock==0
            [pocket] = pocket2find_PL_region(protein,ligand,Rcutoff);
            clear protein_centroid
            protein_centroidA(1,1) = mean(protein(:,1));
            protein_centroidA(1,2) = mean(protein(:,2));
            protein_centroidA(1,3) = mean(protein(:,3));
            
            clear pocket_centroid
            pocket_centroidA(1,1) = mean(pocket(:,1));
            pocket_centroidA(1,2) = mean(pocket(:,2));
            pocket_centroidA(1,3) = mean(pocket(:,3));
            
            
            
            ligand_centroid(1,1)=mean(ligand(:,1));
            ligand_centroid(1,2)=mean(ligand(:,2));
            ligand_centroid(1,3)=mean(ligand(:,3));
            
            lcp_dist_vec = centpolarradi(ligand_centroid,ligand);
            lcpc_dist_vec = dist(ligand_centroid,protein_centroidA);
            
            lcp_dist_max = max(lcp_dist_vec);
            
            %water_mat = watermat(protein,pocket_centroidA);
            protein_surface = [];
            grid_x = 0.1;
            for lyr = 0:0.5:20
                [protein_surface_temp] = pocketsurface(lyr,ligand_centroid,protein_centroidA,lcpc_dist_vec,protein,grid_x);
                
                protein_surface(size(protein_surface,1)+1:size(protein_surface,1)+size(protein_surface_temp,1),:) = protein_surface_temp;
                
            end
            
            protein_surface = unique(protein_surface,'rows');
            protein_surface = sortrows(protein_surface,7);
            
            
            
            clear Grid_Heatmap
            Grid_Heatmap = Heatmap_2a(protein_surface,protein,ligand);
            ligpock=1;
            %end
            %Grid_Heatmap = Heatmap_2a(pocket,protein,ligand);
            
            
            toc
            
            tic
            
            ligand_sig = ligandsig(ligand);
            %indemp1 = 0;
            %indemp2 = 0;
            
            polar_ligand = ligand(ligand(:,4)>=5&ligand(:,4)<=13,:);
            
            
            dist_llc = dist(polar_ligand,ligand_centroid);
            
            
            max_dist_llc = max(dist_llc);
            
            
            %if ligand_sig<3 %|| (indemp1 == 0 && indemp2 == 0)
            pocket_protein_centroid_dist=[];
            pocket_protein_centroid_dist_0=[];
            if isempty(pocket_protein_centroid_dist) && isempty(pocket_protein_centroid_dist_0)
                pocket_centroid = pocket_centroidA;
                pocket_centroid_ave=pocket_centroid;
                pocket_protein_centroid_dist = dist(pocket_centroidA,protein_centroidA);
                pocket_protein_centroid_dist_0 = pocket_protein_centroid_dist;
            end
            lyr_num = 1;
            %lyr_num = 12;
            grid_d = 0.5;
            for lyrx = -20:grid_d:20
                for lyrz = -20:grid_d:20
            %for lyrx = -7.5:grid_d:-7.5
            %    for lyrz = -6:grid_d:-5.5        
                    Panchorx = [lyrx 0 0];
                    Panchorz = [0 0 lyrz];
                    
                    Panchorx = Panchorx + ligand_centroid;
                    Panchorz = Panchorz + ligand_centroid;
                    
                    
                    [pocket_centroid,pocket_centroid_ave] = gridshell_xy(lyrx,lyrz,Panchorx,Panchorz,ligand_centroid,protein_centroidA,protein,protein_surface,grid_d,max_dist_llc);
                    
                    if ~isempty(pocket_centroid)
                        pocket_protein_centroid_dist = dist(pocket_centroid,protein_centroidA);
                        
                    elseif isempty(pocket_centroid)
                        pocket_protein_centroid_dist = pocket_protein_centroid_dist_0;
                    end
                    
                    
                    if ~isempty(pocket_centroid)
                        clear water_mat water_mat_temp
                        water_mat=[];
                        for npc = 1:size(pocket_centroid,1)
                            water_mat_temp = watermat2(protein,ligand_centroid,pocket_centroid(npc,:),lcp_dist_vec);
                            water_mat(size(water_mat,1)+1:size(water_mat,1)+size(water_mat_temp,1),:)=water_mat_temp;
                        end
                        
                        water_mat = unique(water_mat,'rows');
                        water_num=[1:1:size(water_mat,1)]';
                        
                        water_mat(:,7)=water_num;
                        
                        if ~isempty(water_mat) && size(pocket_centroid,1)>0
                            Grid_W_Heatmap = Heatmap_pl_water(water_mat,protein,ligand);
                            [ligand_final2] = DioP_anchor_vf_water(ligand,Grid_W_Heatmap,pocket_centroid,protein,water_mat,protein_centroidA);
                        end
                        
                    end
                    
                    
                    %pocket_lyr_all(lyr_num).lyr = pocket_centroid;
                    %pocket_lyr_all(lyr_num).lyr_ave = pocket_centroid_ave;
                    
                    %if lyr>0
                    %    if dist(pocket_lyr_all(lyr_num).lyr_ave,pocket_centroidA)<=dist(pocket_lyr_all(lyr_num-1).lyr_ave,pocket_centroidA)
                    %        break
                    %    end
                    %end
                    
                    %if abs(lyr - 5.0)<0.1 || abs(lyr - 10.0)<0.1 || abs(lyr - 15.0)<0.1 || abs(lyr - 20.0)<0.1
                    %    [pocket] = pocket2find_PL_region(protein,lig_struc_rec_sel,Rcutoff);
                    %    clear Grid_Heatmap
                    %    Grid_Heatmap = Heatmap_2a(pocket,protein,lig_struc_rec_sel);
                    %end
                    
                    
                    if size(pocket_centroid,1)>0
                        
                        
                        clear ligand_final
                        %[ligand_final] = DioP_anchor_v5(ligand,Grid_Heatmap,pocket_centroid,pocket,protein,protein_centroidA);
                        %[ligand_final] = DioP_anchor_hydpho_LL(ligand,Grid_Heatmap,pocket_centroid,protein,protein_centroidA);
                        [ligand_final1] = DioP_anchor_xy(ligand,Grid_Heatmap,pocket_centroid,protein,protein_centroidA,lcp_dist_max);
                        
                        if ~isempty(water_mat)
                            ligand_final = ligand_final1;
                            for lf2 = 1:size(ligand_final2,2)
                                ligand_final(size(ligand_final1,2)+lf2).dG = ligand_final2(lf2).dG;
                                ligand_final(size(ligand_final1,2)+lf2).structure = ligand_final2(lf2).structure;
                                ligand_final(size(ligand_final1,2)+lf2).repulse = ligand_final2(lf2).repulse;
                                ligand_final(size(ligand_final1,2)+lf2).procentroid_dist = ligand_final2(lf2).procentroid_dist;
                            end
                        else
                            ligand_final = ligand_final1;
                        end
                        
                        toc
                        
                        clear ligand_final_score
                        ligand_final_score=zeros(size(ligand_final,2),3);
                        lfn=1;
                        for lf=1:size(ligand_final,2)
                            if lyr_num == 1
                                
                                ligand_final_score(lfn,1) = ligand_final(lf).dG;
                                ligand_final_score(lfn,2) = lf;
                                ligand_final_score(lfn,3) = ligand_final(lf).procentroid_dist;
                                lfn=lfn+1;
                            elseif lyr_num > 1 && ~isempty(ligand_docked_dist_standard)
                                if ligand_final(lf).procentroid_dist >= ligand_docked_dist_standard
                                    ligand_final_score(lfn,1) = ligand_final(lf).dG;
                                    ligand_final_score(lfn,2) = lf;
                                    ligand_final_score(lfn,3) = ligand_final(lf).procentroid_dist;
                                    lfn=lfn+1;
                                elseif ligand_final(lf).procentroid_dist < ligand_docked_dist_standard
                                    ligand_final_score(lfn,1) = 100;
                                    ligand_final_score(lfn,2) = lf;
                                    ligand_final_score(lfn,3) = ligand_final(lf).procentroid_dist;
                                    lfn=lfn+1;
                                    
                                end
                            elseif lyr_num > 1 && isempty(ligand_docked_dist_standard)
                                ligand_final_score(lfn,1) = 100;
                                ligand_final_score(lfn,2) = lf;
                                ligand_final_score(lfn,3) = ligand_final(lf).procentroid_dist;
                                lfn=lfn+1;
                            end
                        end
                        
                        ligand_docked_dist_standard= max(ligand_final_score(:,3));
                        
                        ligand_final_score(ligand_final_score(:,1)==0,:)=[];
                        
                        ligand_final_score = sortrows(ligand_final_score,[1]);
                        clear ligand_select_X1
                        ligand_select_X1=struct([]);
                        for lsxn=1:size(ligand_final,2)
                            
                            ligand_select_X1(lsxn).coordinate = ligand_final(ligand_final_score(lsxn,2)).structure;
                            ligand_select_X1(lsxn).dG_X1 = ligand_final(ligand_final_score(lsxn,2)).dG;
                            ligand_select_X1(lsxn).repulse = ligand_final(ligand_final_score(lsxn,2)).repulse;
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        clear dG_select_X1 dG_select_X2 dG_select_pass_X1 dG_select_pass_X2 Final_result_X1 Final_result_X2
                        dG_select_X1=[];
                        dG_select_pass_X1=[];
                        Final_result_X1=[];
                        for i = 1:size(ligand_select_X1,2)
                            if ~isempty(ligand_select_X1(i).dG_X1)
                                if ligand_select_X1(i).repulse<=60
                                    dG_select_X1(i,1)=ligand_select_X1(i).dG_X1;
                                    dG_select_X1(i,2)=ligand_select_X1(i).repulse;
                                elseif ligand_select_X1(i).repulse>60
                                    dG_select_X1(i,1)=0;
                                    dG_select_X1(i,2)=0;
                                end
                            elseif isempty(ligand_select_X1(i).dG_X1)
                                dG_select_X1(i,1)=0;
                                dG_select_X1(i,2)=0;
                            end
                        end
                        
                        if size(dG_select_X1,1)>0
                            dG_select_X1(:,3)=1:size(dG_select_X1,1);
                            dG_select_X1(dG_select_X1(:,1)==0,:)=[];
                            
                            dG_select_X1 = sortrows(dG_select_X1,[1]);
                            
                            if size(dG_select_X1,1)>50
                                dG_select_pass_X1 = dG_select_X1(1:50,:);
                            else
                                dG_select_pass_X1 = dG_select_X1;%(dG_select_X1(:,1)<0,:);
                            end
                            
                            dG_select_pass_X1(dG_select_pass_X1(:,1)==0,:)=[];
                            
                            dG_select_pass_X1 = sortrows(dG_select_pass_X1,[1]);
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            mkdir(strcat('F:\MTdockPMF_water_Case_xyF1\Trypsin_Benzamidine\',ligand_name,'\',num2str(lyrx),'_',num2str(lyrz),'_Angstrom_layer','\output\'));
                            clear lig_struc_record lig_struc_rec_sel
                            for i = 1:size(dG_select_pass_X1,1)
                                %[ligand_final_str_X1,ligand_final_dG_X1,ligand_final_repulse_num_X1] = Erelax_rotaC1v2(ligand_select_X1,dG_select_pass_X1,i,pocket,vdw_d);
                                
                                [ligand_final_str_X1,ligand_final_dG_X1,ligand_final_repulse_num_X1] = Erelax_NRG_water(ligand_select_X1,dG_select_pass_X1,i,protein,water_mat,vdw_d);
                                RMSD_final = RMSD_cal(ligand_final_str_X1,ligand_origin);
                                ligand_str_X1_centroid(1,1) = mean(ligand_final_str_X1(:,1));
                                ligand_str_X1_centroid(1,2) = mean(ligand_final_str_X1(:,2));
                                ligand_str_X1_centroid(1,3) = mean(ligand_final_str_X1(:,3));
                                
                                Final_result_X1(i,1)=ligand_final_dG_X1(1,1);
                                Final_result_X1(i,2)=ligand_final_repulse_num_X1;
                                Final_result_X1(i,3)=RMSD_final;
                                Final_result_X1(i,4)=0;
                                Final_result_X1(i,5)=abs(ligand_str_X1_centroid(1,1)-ligand_centroid(1,1));
                                Final_result_X1(i,6)=abs(ligand_str_X1_centroid(1,3)-ligand_centroid(1,3));
                                %if abs(lyr - 4.5)<0.1 || abs(lyr - 9.5)<0.1 || abs(lyr - 14.5)<0.1 || abs(lyr - 19.5)<0.1
                                %    lig_struc_record(i).input_struc = ligand_final_str_X1;
                                %end
                                
                                if Final_result_X1(i,2)<=8 && abs(Final_result_X1(i,5)-abs(lyrx)) <=2*grid_d && abs(Final_result_X1(i,6)-abs(lyrz)) <=2*grid_d
                                    Final_result_X1(i,4)=1;
                                    input_struc = ligand_final_str_X1;
                                    
                                    
                                    listC = strcat(ligand_name,'_X1');
                                    output_dir = strcat('F:\MTdockPMF_water_Case_xyF1\Trypsin_Benzamidine\',ligand_name,'\',num2str(lyrx),'_',num2str(lyrz),'_Angstrom_layer','\output\');
                                    
                                    out_X1 = lig_layer_xy_output(input_struc,listC,lyrx,lyrz,i,output_dir);
                                end
                                
                            end
                            
                            
                            %if abs(lyr - 4.5)<0.1 || abs(lyr - 9.5)<0.1 || abs(lyr - 14.5)<0.1 || abs(lyr - 19.5)<0.1
                            %    lyr_measure = ((lyr-grid_d)+(lyr+grid_d))/2;
                            %    clear lig_struc_record_temp
                            %    ti=1;
                            %    for i = 1:size(lig_struc_record,2)
                            %        if ~isempty(lig_struc_record(i).input_struc)
                            %            lig_struc_record_temp(ti).input_struc = lig_struc_record(i).input_struc;
                            %            ti=ti+1;
                            %        end
                            %    end
                            %    clear lig_struc_record
                            %    lig_struc_record = lig_struc_record_temp;
                            %    dG_max_ind = find(abs(Final_result_X1(:,5)-lyr_measure) == min(abs(Final_result_X1(:,5)-lyr_measure)));
                            %    lig_struc_rec_sel = lig_struc_record(dG_max_ind).input_struc;
                            %    clear lig_struc_record
                            %end
                            
                            
                            Final_result_X1(Final_result_X1(:,4)==0,:)=[];
                            Final_result_X1(:,4)=[];
                            
                            
                            if size(Final_result_X1,1) >0
                                Final_result=[Final_result_X1
                                    ];
                                
                                data_layer_xy_output(Final_result,lyrx,lyrz,output_dir,code_dir);
                                lyr_num=lyr_num+1
                                %cd('C:\Users\John Zheng\Documents\MATLAB\Heatmap_Docking');
                            elseif size(Final_result_X1,1) ==0
                                lyr_num=lyr_num+1
                            end
                            %Final_result_data(jpz).data=Final_result;
                        elseif size(dG_select_X1,1) == 0 && lyr_num>1
                            lyr_num=lyr_num+1
                        elseif size(dG_select_X1,1) == 0 && lyr_num==1
                            continue
                            
                        end
                        
                        
                        %end
                        
                    elseif size(pocket_centroid,1) == 0
                        continue
                    end
                    
                end
            end
        end
    end
end
