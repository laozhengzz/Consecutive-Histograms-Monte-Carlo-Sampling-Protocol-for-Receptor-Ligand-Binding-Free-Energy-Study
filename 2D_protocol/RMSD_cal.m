function rmsd = RMSD_cal(data1,data2)

for i =1:size(data1,1)
    vari(i,1)=((data1(i,1)-data2(i,1)).^2+(data1(i,2)-data2(i,2)).^2+(data1(i,3)-data2(i,3)).^2);
end
rmsd=sqrt(sum(vari)/size(vari,1));
    
    