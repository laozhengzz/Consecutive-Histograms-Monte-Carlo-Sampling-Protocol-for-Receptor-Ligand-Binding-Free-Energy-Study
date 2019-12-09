clear all

home_dir = strcat('C:\Users\John Zheng\Documents\pdb_test\');
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
            [protein,protein_Namelist] = protein_define_final(protein_dir);
        end
    end
end

%39
for zh=1:size_list1(1,1)
%for zh=7:7
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
            %ligand(ligand(:,4)==6,4)=5;
            
            %boundary_id = [ strcat(ligand_name(1:end-5),'-blob')];
            %boundary_name = [ strcat(boundary_id,'.txt')];
            %boundary_dir = [ strcat(home_dir,pdbid,'\',boundary_name)];
            
            %[region]= boundary_define_simp(boundary_dir);
            
            %if ligpock==0
            [pocket] = pocket2find_PL_region(protein,ligand,Rcutoff);
            clear Grid_Heatmap
            %Grid_Heatmap = Heatmap_2a(pocket,protein,ligand);
            ligpock=1;
            %end
            %Grid_Heatmap = Heatmap_2a(pocket,protein,ligand);
            
            
            toc
            
            
            
            tic
            
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
            
            ligand_sig = ligandsig(ligand);
            %indemp1 = 0;
            %indemp2 = 0;
            
            %if ligand_sig<3 %|| (indemp1 == 0 && indemp2 == 0)
            pocket_centroid_all=[];
            
            lyr_num = 1;
            for lyr = 0:0.4:2.8
                if lyr == 0
                    pocket_centroid = pocket_centroidA;
                    pocket_protein_centroid_dist = dist(pocket_centroidA,protein_centroidA);
                    pocket_protein_centroid_dist_0 = pocket_protein_centroid_dist;
                elseif lyr > 0
                    pocket_centroid = gridshell(lyr,pocket_centroidA,protein_centroidA,pocket_protein_centroid_dist,protein);
                    if ~isempty(pocket_centroid)
                        pocket_protein_centroid_dist = dist(pocket_centroid,protein_centroidA);
                    elseif isempty(pocket_centroid)
                        pocket_protein_centroid_dist = pocket_protein_centroid_dist_0;
                    end
                end
                
                pocket_centroid_all(size(pocket_centroid_all,1)+1:size(pocket_centroid_all,1)+size(pocket_centroid,1),:) = pocket_centroid;
                
                
            end
            
        end
    end
end
