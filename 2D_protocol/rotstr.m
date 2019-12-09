function structure_new = rotstr(structure,axis,degree)
structure_new = zeros(size(structure,1),10);
%structure_sph(:,1) = sqrt(structure(:,1).^2 + structure(:,2).^2 + structure(:,3).^2);
%structure_sph(:,2) = acos(structure(:,3)./sqrt(structure(:,1).^2 + structure(:,2).^2 + structure(:,3).^2));

if axis == 1
    p1=1;
    p2=2;
    p3=3;
elseif axis == 2
    p1=2;
    p2=1;
    p3=3;
elseif axis == 3
    p1=3;
    p2=1;
    p3=2;
end
    



for i = 1:size(structure,1)
    structure_sph(i,1)=sqrt(structure(i,p2).^2 + structure(i,p3).^2);
    if structure(i,p2)>0
        structure_sph(i,2)=atan(structure(i,p3)/structure(i,p2));
    elseif structure(i,p2)<0
        structure_sph(i,2)=pi+atan(structure(i,p3)/structure(i,p2));
    elseif structure(i,p2)==0 && structure(i,p3)>0
        structure_sph(i,2)=pi/2;
    elseif structure(i,p2)==0 && structure(i,p3)<0
        structure_sph(i,2)=-pi/2;
    end
    
    structure_sph(i,2)=structure_sph(i,2)+degree;
    
    structure_new(i,p1) = structure(i,p1);
    structure_new(i,p2) = structure_sph(i,1)*cos(structure_sph(i,2));
    structure_new(i,p3) = structure_sph(i,1)*sin(structure_sph(i,2));
    
end

structure_new(:,4:10)=structure(:,4:10);