function pointgroupf = findPCB1_correct(lyr,grid_d,centroid,protein_centroidA,pocket,protein_structure)
%centroid = pocket_centroidA;
%protein_structure = protein;
vec_pro_poc = centroid(1,1:3) - protein_centroidA(1,1:3);
vec_poc_pro = protein_centroidA(1,1:3) - centroid(1,1:3);
protein_centroidAdj = vec_poc_pro;

rpca = sqrt(protein_centroidAdj(1,1)^2+protein_centroidAdj(1,2)^2+protein_centroidAdj(1,3)^2);
%
thetapca = acos(protein_centroidAdj(1,3)/rpca);
if protein_centroidAdj(1,2)>0
    phipca = acos(protein_centroidAdj(1,1)/(rpca*sin(thetapca)));
elseif protein_centroidAdj(1,2)<0
    phipca = -acos(protein_centroidAdj(1,1)/(rpca*sin(thetapca)));
end
clear sphpointgroup carpointgroup vec_carp_poc
%sphpointgroup(1,:) = [rpca thetapca phipca];
rpca = 5;
scount = 1;
for i = 1:10
    for j = 1:10
        for k = 1:4
            if k == 1
                sphpointgroup(scount,1)=rpca;
                sphpointgroup(scount,2)=thetapca + (i-1)*(18*pi/180);
                sphpointgroup(scount,3)=phipca + (j-1)*(18*pi/180);
                scount = scount+1;
            elseif k == 2
                sphpointgroup(scount,1)=rpca;
                sphpointgroup(scount,2)=thetapca - (i-1)*(18*pi/180);
                sphpointgroup(scount,3)=phipca - (j-1)*(18*pi/180);
                scount = scount+1;
            elseif k == 3
                sphpointgroup(scount,1)=rpca;
                sphpointgroup(scount,2)=thetapca + (i-1)*(18*pi/180);
                sphpointgroup(scount,3)=phipca - (j-1)*(18*pi/180);
                scount = scount+1;
            elseif k == 4
                sphpointgroup(scount,1)=rpca;
                sphpointgroup(scount,2)=thetapca - (i-1)*(18*pi/180);
                sphpointgroup(scount,3)=phipca + (j-1)*(18*pi/180);
                scount = scount+1;    
                
                
            end
        end
    end
end

sphpointgroup = unique(sphpointgroup,'rows');

for i = 1:size(sphpointgroup,1)
    carpointgroup(i,1)=sphpointgroup(i,1)*sin(sphpointgroup(i,2))*cos(sphpointgroup(i,3)) + centroid(1,1);
    carpointgroup(i,2)=sphpointgroup(i,1)*sin(sphpointgroup(i,2))*sin(sphpointgroup(i,3)) + centroid(1,2);
    carpointgroup(i,3)=sphpointgroup(i,1)*cos(sphpointgroup(i,2)) + centroid(1,3);
end


%scatter3(carpointgroup(:,1),carpointgroup(:,2),carpointgroup(:,3))
%axis equal

vec_carp_poc(:,1) = centroid(1,1) - carpointgroup(:,1);
vec_carp_poc(:,2) = centroid(1,2) - carpointgroup(:,2);
vec_carp_poc(:,3) = centroid(1,3) - carpointgroup(:,3);

RCvec_carp_pok = zeros(size(carpointgroup,1),2);


    
    
for i = 1:size(protein_structure,1)
    
    vec_carp_pok(:,1) = protein_structure(i,1) - carpointgroup(:,1);
    vec_carp_pok(:,2) = protein_structure(i,2) - carpointgroup(:,2);
    vec_carp_pok(:,3) = protein_structure(i,3) - carpointgroup(:,3);
    
    for j = 1:size(carpointgroup,1)
        if dist(carpointgroup(j,1:3),protein_structure(i,1:3))<2.6
            RCvec_carp_pok(j,2)=1;
            continue
        end
            
            
        cosangle = dot(vec_carp_poc(j,:),vec_carp_pok(j,:))/(norm(vec_carp_poc(j,:))*norm(vec_carp_pok(j,:)));
        %if cosangle<0
        %    continue
        %end
        Rvec_carp_pok = norm(vec_carp_pok(j,:));
        sinangle = sqrt(1-cosangle^2);
        Rvec_carp_pok_s = Rvec_carp_pok*sinangle;
        Rvec_carp_pok_sdot(j,1) = Rvec_carp_pok_s;
        if Rvec_carp_pok_s<2.6
            RCvec_carp_pok(j,1)=RCvec_carp_pok(j,1)+1;
        end
    end
    
    
end

[pm1,pm2]=find(RCvec_carp_pok(:,1) <= min(RCvec_carp_pok(RCvec_carp_pok(:,2)~=1,1))+0);
%[pm1,pm2]=find(RCvec_carp_pok(:,1) <= min(RCvec_carp_pok(:,1))+5);

clear innerpoint1 vec_carp_poc1 vec_carp_pok1
RCvec_carp_pok1 = zeros(size(pm1,1),1);
if size(pm1,1)>1
    
    for ia = 1:size(pm1,1)
        innerpoint1(ia,:) = carpointgroup(pm1(ia,1),:);
    end
    
    vec_carp_poc1(:,1) = centroid(1,1) - innerpoint1(:,1);
    vec_carp_poc1(:,2) = centroid(1,2) - innerpoint1(:,2);
    vec_carp_poc1(:,3) = centroid(1,3) - innerpoint1(:,3);

    for i = 1:size(pocket,1)
        
        vec_carp_pok1(:,1) = pocket(i,1) - innerpoint1(:,1);
        vec_carp_pok1(:,2) = pocket(i,2) - innerpoint1(:,2);
        vec_carp_pok1(:,3) = pocket(i,3) - innerpoint1(:,3);
        
        for j = 1:size(innerpoint1,1)
            cosangle = dot(vec_carp_poc1(j,:),vec_carp_pok1(j,:))/(norm(vec_carp_poc1(j,:))*norm(vec_carp_pok1(j,:)));
            cosangledot(j,1) = cosangle;
            Rvec_carp_pok1 = norm(vec_carp_pok1(j,:));
            sinangle = sqrt(1-cosangle^2);
            Rvec_carp_pok_s = Rvec_carp_pok1*sinangle;
            %if Rvec_carp_pok_s<2.6
                RCvec_carp_pok1(j,1)=RCvec_carp_pok1(j,1)+Rvec_carp_pok_s^2;
            %end
        end
        
        
        
    end
end

[pm11,pm22]=find(RCvec_carp_pok1 > max(RCvec_carp_pok1)-0.01);
    
nia=1;
for ia = 1:size(pm11,1)
    innerpoint = carpointgroup(pm1(pm11(ia,1),1),:);
    vec1 = centroid(1,:) - innerpoint(1,:);
    vec2 = -centroid(1,:) + innerpoint(1,:);
    
    ancpoint1 = centroid(1,:) + vec1;
    ancpoint2 = centroid(1,:) + vec2;
    
    if min(dist(protein_structure(:,1:3),ancpoint1(1,1:3)))<min(dist(protein_structure(:,1:3),ancpoint2(1,1:3)))
        centroidAdj = vec2;
    elseif min(dist(protein_structure(:,1:3),ancpoint1(1,1:3)))>min(dist(protein_structure(:,1:3),ancpoint2(1,1:3)))
        centroidAdj = vec1;
    end
    
    %centroidAdj = -centroid(1,:) + innerpoint(1,:);
    
    rpca = sqrt(centroidAdj(1,1)^2+centroidAdj(1,2)^2+centroidAdj(1,3)^2);
    thetapca = acos(centroidAdj(1,3)/rpca);
    if centroidAdj(1,2)>0
        phipca = acos(centroidAdj(1,1)/(rpca*sin(thetapca)));
    elseif centroidAdj(1,2)<0
        phipca = -acos(centroidAdj(1,1)/(rpca*sin(thetapca)));
    end
    
    anglechg= acos((2*rpca^2-(5*grid_d/lyr)^2)/(2*rpca^2));
    
    clear sphpointgroupf carpointgroupf
    %sphpointgroup(1,:) = [rpca thetapca phipca];
    scount = 1;
    for i = 1:5
        for j = 1:5
            for k = 1:4
                if k == 1
                    sphpointgroupf(scount,1)=lyr;
                    sphpointgroupf(scount,2)=thetapca + (i-1)*anglechg;
                    sphpointgroupf(scount,3)=phipca + (j-1)*anglechg;
                    scount = scount+1;
                elseif k == 2
                    sphpointgroupf(scount,1)=lyr;
                    sphpointgroupf(scount,2)=thetapca - (i-1)*anglechg;
                    sphpointgroupf(scount,3)=phipca - (j-1)*anglechg;
                    scount = scount+1;
                elseif k == 3
                    sphpointgroupf(scount,1)=lyr;
                    sphpointgroupf(scount,2)=thetapca + (i-1)*anglechg;
                    sphpointgroupf(scount,3)=phipca - (j-1)*anglechg;
                    scount = scount+1;
                elseif k == 4
                    sphpointgroupf(scount,1)=lyr;
                    sphpointgroupf(scount,2)=thetapca - (i-1)*anglechg;
                    sphpointgroupf(scount,3)=phipca + (j-1)*anglechg;
                    scount = scount+1;
                    
                end
            end
        end
    end
    
    sphpointgroupf = unique(sphpointgroupf,'rows');
    
    for i = 1:size(sphpointgroupf,1)
        carpointgroupf(i,1)=sphpointgroupf(i,1)*sin(sphpointgroupf(i,2))*cos(sphpointgroupf(i,3)) + centroid(1,1);
        carpointgroupf(i,2)=sphpointgroupf(i,1)*sin(sphpointgroupf(i,2))*sin(sphpointgroupf(i,3)) + centroid(1,2);
        carpointgroupf(i,3)=sphpointgroupf(i,1)*cos(sphpointgroupf(i,2)) + centroid(1,3);
    end
    
    
if size(protein_structure,1)>0
    for jz=1:1:size(protein_structure,1)
        carpointgroupf(sqrt( (carpointgroupf(:,1)-protein_structure(jz,1)).^2+(carpointgroupf(:,2)-protein_structure(jz,2)).^2+(carpointgroupf(:,3)-protein_structure(jz,3)).^2 )< 2.6,:)=[];
    end
end

    
    
    pointgroupf(ia).carpointgroupf = carpointgroupf;
    

%scatter3(carpointgroupf(:,1),carpointgroupf(:,2),carpointgroupf(:,3))
%axis equal

end
    
    
    




