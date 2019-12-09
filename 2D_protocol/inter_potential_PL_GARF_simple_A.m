%function [ZE_part_PL,ZE_part_PL2,ZN_part_PL,ZN_part_PL2] = inter_potential(protein_starterA,protein_starterB,vdw_d)
function [Part_Matr_com,CNN1,VN1] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d)
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2;
const1=halv1*2;
column_num=100;
Part_Matr_com=0;
Part_Matr_comA=0;
for i=1:1:size(protein_starterA,1)
    
    a = protein_starterA(i,4);
    b = protein_starterB(i,4);
    
    dist=norm(protein_starterB(i,1:3)-protein_starterA(i,1:3)); %%cf_vdw_new=sum(cf_vdw);
    sd=[];
    sd=find(abs(vdw_d(:,1)-dist)<=rd/2);   %%cf_vdw_new=1/cf_vdw_new;
    
    if isempty(sd)
        continue
    end
    
    order_num=(a-1)*34+b+1;
    %c3n2
    if (sd>halva) & (sd<length(vdw_d)-halva) %end
        com_vdw=vdw_d(sd-halva+1:sd+halva,order_num);
        %cf_vdw=repmat(vdw_f(sd-halva+1:sd+halva,order_num),round(800/(2*halva)),1);
    elseif sd>=length(vdw_d)-halva
        com_vdw=vdw_d(sd-(2*halva-1):sd,order_num);
        %cf_vdw=repmat(vdw_f(sd-(2*halva-1):sd,order_num),round(800/(2*halva)),1);
    elseif sd<=halva
        com_vdw=vdw_d(sd:sd+(2*halva-1),order_num);
        %cf_vdw=repmat(vdw_f(sd:sd+(2*halva-1),order_num),round(800/(2*halva)),1);
    end
    
    com_vdw_new=sum(com_vdw);
    
    if com_vdw_new>0
        Part_Matr_com=Part_Matr_com+log(com_vdw_new);
        
        %Eij(i,1) = -0.5918*log(com_vdw_new);
    end
    
    
    
end




%Part_Matr_com = Part_Matr_com*1.4;



CNN1=size(protein_starterA,1);

VN1 = CNN1*log(const);














