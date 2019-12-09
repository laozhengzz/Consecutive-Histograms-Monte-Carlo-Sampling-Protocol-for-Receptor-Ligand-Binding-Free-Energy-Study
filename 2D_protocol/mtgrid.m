function Emt_point = mtgrid(pocket_structure,protein_whole,grid_point,grid_type,Paraset,grid_d,Dist_Mat_prep,jz)


b = grid_type;
rd = 0.005;
hnstate = grid_d/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2+1;
const1=halv1*2+1;
column_num=100;
Part_Matr_com=0;
Part_Matr_comA=0;

Emt_point = zeros(size(grid_point,1),4);
%Dist_Mat = zeros(size(grid_point,1),const);

%Dist_Mat_prep = zeros(size(grid_point,1),const);

%Dist_Mat_prep(:,halva+1)=0;
%for j = 1:halva
%    Dist_Mat_prep(:,halva+1-j)=0-j*rd;
%    Dist_Mat_prep(:,halva+1+j)=j*rd;
%end

input_structure = protein_whole(sqrt((protein_whole(:,1)-pocket_structure(jz,1)).^2+(protein_whole(:,2)-pocket_structure(jz,2)).^2+(protein_whole(:,3)-pocket_structure(jz,3)).^2)>=0 & sqrt((protein_whole(:,1)-pocket_structure(jz,1)).^2+(protein_whole(:,2)-pocket_structure(jz,2)).^2+(protein_whole(:,3)-pocket_structure(jz,3)).^2)<6,:);
CNN1=size(input_structure,1);
VN1 = CNN1*log(const);
%for h = 1:size(grid_point,1)
%    pocket_structure = input_structure(sqrt((grid_point(h,1)-input_structure(:,1)).^2+(grid_point(h,2)-input_structure(:,2)).^2+(grid_point(h,3)-input_structure(:,3)).^2)<5,:);

Part_Matr_com=zeros(size(grid_point,1),1);
for i=1:1:size(input_structure,1)
    a = input_structure(i,4);
    order_num=(a-1)*34+b+1;
    dist=sqrt((grid_point(:,1)-input_structure(i,1)).^2+(grid_point(:,2)-input_structure(i,2)).^2+(grid_point(:,3)-input_structure(i,3)).^2); %%cf_vdw_new=sum(cf_vdw);
    
    Dist_Mat = repmat(dist,1,const);
    
    Dist_Mat = Dist_Mat+Dist_Mat_prep;
    
    %
    com_vdw = exp((Paraset(order_num,1).*(Paraset(order_num,2)./Dist_Mat).^Paraset(order_num,3)+Paraset(order_num,4).*(Paraset(order_num,5)./Dist_Mat).^Paraset(order_num,6))./-0.5918);
    %
    com_vdw_new=log(sum(com_vdw')');
    Part_Matr_com = Part_Matr_com + com_vdw_new;
    
    
end



Emt_point(:,1:3)=grid_point(:,1:3);
Emt_point(:,4) = -0.5918.*(Part_Matr_com - VN1);

Emt_point = sortrows(Emt_point,[4,1,2,3]);

EPosi_num = size(Emt_point(Emt_point(:,4)>0,:),1);

onethird = floor(size(Emt_point,1)/3);

if Emt_point(1,4)>0 || EPosi_num/size(Emt_point,1)<1/3
    Emt_point = Emt_point(1:onethird,:);
elseif EPosi_num/size(Emt_point,1)>=1/3
    Emt_point = Emt_point(Emt_point(:,4)<0,:);
end


%end





