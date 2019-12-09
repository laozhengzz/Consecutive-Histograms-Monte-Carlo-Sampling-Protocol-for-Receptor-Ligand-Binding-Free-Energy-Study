%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = lig_layer_output(input_struc,list,lyr,asdf,output_dir)

nonpolar_ligand_n_n=input_struc;

size_nonpolar_ligand_n_n=size(nonpolar_ligand_n_n);
%size_nonpolar_protein_n_n=size(nonpolar_protein_n_n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ligand


for i=1:1:size_nonpolar_ligand_n_n(1,1)
    if ((nonpolar_ligand_n_n(i,4)>=1)&&(nonpolar_ligand_n_n(i,4)<=4))||((nonpolar_ligand_n_n(i,4)>=26)&&(nonpolar_ligand_n_n(i,4)<=34))||(nonpolar_ligand_n_n(i,4)==19)
        nonpolar_ligand_n_n_char(i,:)='C ';
    elseif ((nonpolar_ligand_n_n(i,4)>=5)&&(nonpolar_ligand_n_n(i,4)<=9))
        nonpolar_ligand_n_n_char(i,:)='O ';
    elseif ((nonpolar_ligand_n_n(i,4)>=10)&(nonpolar_ligand_n_n(i,4)<=13))
        nonpolar_ligand_n_n_char(i,:)='N ';
    elseif (nonpolar_ligand_n_n(i,4)>=20)&(nonpolar_ligand_n_n(i,4)<=21)
        nonpolar_ligand_n_n_char(i,:)='S ';
    elseif (nonpolar_ligand_n_n(i,4)==15)
        nonpolar_ligand_n_n_char(i,:)='F ';
    elseif (nonpolar_ligand_n_n(i,4)==16)
        nonpolar_ligand_n_n_char(i,:)='Cl';
    elseif (nonpolar_ligand_n_n(i,4)==17)
        nonpolar_ligand_n_n_char(i,:)='Br';
    elseif (nonpolar_ligand_n_n(i,4)==18)
        nonpolar_ligand_n_n_char(i,:)='I ';
    elseif (nonpolar_ligand_n_n(i,4)==14)
        nonpolar_ligand_n_n_char(i,:)='P ';
    elseif (nonpolar_ligand_n_n(i,4)>=22)&(nonpolar_ligand_n_n(i,4)<=25)
        nonpolar_ligand_n_n_char(i,:)='H ';
    else
        nonpolar_ligand_n_n_char(i,:)='? ';
    end
end


k=1;
l=1;


%fid=fopen('np_l_n_na.txt','wt')
%for a=1:1:size_nonpolar_ligand_n_n_charge(1,1)
%    fprintf(fid,'% 5.2d', nonpolar_ligand_n_n_num(a,1));
%    fprintf(fid,'\npair.xyz% 5.2d\n', nonpolar_ligand_n_n_charge(a,1));
%    l=a
%    for i=1:1:size_nonpolar_ligand_n_n(1,1)
%        if nonpolar_ligand_n_n(i,10)==nonpolar_ligand_n_n(l,10)
%            fprintf(fid,'% 5.2c', nonpolar_ligand_n_n_char(i,1));
%            fprintf(fid,'% 5d% 5d% 5d\n', nonpolar_ligand_n_n(i,1),nonpolar_ligand_n_n(i,2),nonpolar_ligand_n_n(i,3));
%        elseif nonpolar_ligand_n_n(i,10)~=nonpolar_ligand_n_n(l,10)
%            fprintf(fid,'=== === ===\n');
%            l=i
%        end
%    end
%end

clear fid
%mkdir(strcat('E:\OSDP1\new\',list));
mkdir(strcat(output_dir));
%script1=[strcat('E:\OSDP\all\',list,'\',list,'_',num2str(asdf),'.txt')];

fid=fopen(strcat(output_dir,'\',list,'_',num2str(lyr),'_',num2str(asdf),'.xyz'),'wt');



%a=1
%fprintf(fid,'% 5.2d', list);
fprintf(fid,'% 5.2d', size_nonpolar_ligand_n_n(1,1));
fprintf(fid,'\npair.xyz% 5.2d\n', 1);
for i=1:1:size_nonpolar_ligand_n_n(1,1)
    
    
    fprintf(fid,'% 5s', nonpolar_ligand_n_n_char(i,:));
    fprintf(fid,'  % 5d  % 5d  % 5d\n', nonpolar_ligand_n_n(i,1),nonpolar_ligand_n_n(i,2),nonpolar_ligand_n_n(i,3));
    
end


fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out=fid;









