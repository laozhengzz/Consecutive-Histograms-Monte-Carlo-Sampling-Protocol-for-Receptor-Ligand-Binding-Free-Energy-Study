%clear all
function protein_list_num = protein_list_final_AF_format(protein_Namelist)


%load torsion_list_final_AF
for i = 1:size(protein_Namelist,1)
    torsion_name1 = protein_Namelist{i,1};
    clear torsion_name_num1
    if length(torsion_name1)<=3
        length_torsion_name1 = length(torsion_name1);
    elseif length(torsion_name1)>3
        length_torsion_name1 = 3;
    end
    for j = 1:length_torsion_name1
        if j == 1 && length(torsion_name1) >=3
            if strcmp(torsion_name1(j:j+2),'OXT')
                torsion_name_num1(j) = '9';
                break
            end
        end
        if strcmp(torsion_name1(j),'C') ==1
            torsion_name_num1(j) = '1';
        elseif strcmp(torsion_name1(j),'N')==1
            torsion_name_num1(j) = '2';
        elseif strcmp(torsion_name1(j),'O')==1
            torsion_name_num1(j) = '3';
        elseif strcmp(torsion_name1(j),'S')==1
            torsion_name_num1(j) = '4';
        elseif strcmp(torsion_name1(j),'A')==1
            torsion_name_num1(j) = '1';
        elseif strcmp(torsion_name1(j),'B')==1
            torsion_name_num1(j) = '2';
        elseif strcmp(torsion_name1(j),'D')==1
            torsion_name_num1(j) = '3';    
        elseif strcmp(torsion_name1(j),'E')==1
            torsion_name_num1(j) = '4';    
        elseif strcmp(torsion_name1(j),'G')==1
            torsion_name_num1(j) = '5';        
        elseif strcmp(torsion_name1(j),'H')==1
            torsion_name_num1(j) = '6';            
        elseif strcmp(torsion_name1(j),'Z')==1
            torsion_name_num1(j) = '7';    
        elseif strcmp(torsion_name1(j),'1')==1
            torsion_name_num1(j) = '1';    
        elseif strcmp(torsion_name1(j),'2')==1
            torsion_name_num1(j) = '2';
        elseif strcmp(torsion_name1(j),'3')==1
            torsion_name_num1(j) = '3';
        else 
            torsion_name_num1(j) = '0';
            
            
        end
    end
    torsion_name1 = protein_Namelist{i,2};
    
    if strcmp(torsion_name1, 'ALA')==1
        torsion_res_num1 = 101;
    elseif strcmp(torsion_name1, 'ARG')==1
        torsion_res_num1 = 102;
    elseif strcmp(torsion_name1, 'ASN')==1 || strcmp(torsion_name1, 'ASX')==1
        torsion_res_num1 = 103;
    elseif strcmp(torsion_name1, 'ASP')==1
        torsion_res_num1 = 104;
    elseif strcmp(torsion_name1, 'CYS')==1
        torsion_res_num1 = 105;
    elseif strcmp(torsion_name1, 'GLN')==1
        torsion_res_num1 = 106;
    elseif strcmp(torsion_name1, 'GLU')==1 || strcmp(torsion_name1, 'GLX')==1
        torsion_res_num1 = 107;
    elseif strcmp(torsion_name1, 'GLY')==1
        torsion_res_num1 = 108;
    elseif strcmp(torsion_name1, 'HIS')==1 || strcmp(torsion_name1, 'HID')==1 || strcmp(torsion_name1, 'HIE')==1 || strcmp(torsion_name1, 'HIP')==1
        torsion_res_num1 = 109;
    elseif strcmp(torsion_name1, 'ILE')==1
        torsion_res_num1 = 110;
    elseif strcmp(torsion_name1, 'LEU')==1
        torsion_res_num1 = 111;
    elseif strcmp(torsion_name1, 'LYS')==1
        torsion_res_num1 = 112;
    elseif strcmp(torsion_name1, 'MET')==1
        torsion_res_num1 = 113;
    elseif strcmp(torsion_name1, 'PHE')==1
        torsion_res_num1 = 114;
    elseif strcmp(torsion_name1, 'PRO')==1
        torsion_res_num1 = 115;
    elseif strcmp(torsion_name1, 'SER')==1
        torsion_res_num1 = 116;
    elseif strcmp(torsion_name1, 'THR')==1
        torsion_res_num1 = 117;
    elseif strcmp(torsion_name1, 'TRP')==1
        torsion_res_num1 = 118;
    elseif strcmp(torsion_name1, 'TYR')==1
        torsion_res_num1 = 119;
    elseif strcmp(torsion_name1, 'VAL')==1
        torsion_res_num1 = 120;
    else
        torsion_res_num1 = 121;
    end
    

    
    
    
    
    protein_list_num(i,1) = str2num(torsion_name_num1);
    protein_list_num(i,2) = torsion_res_num1;
    
end
            
            
            
            
            
            
    