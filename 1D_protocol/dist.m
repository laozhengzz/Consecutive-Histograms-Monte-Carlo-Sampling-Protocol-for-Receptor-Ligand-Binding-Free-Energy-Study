function dist1 = dist(point1,point2)

if size(point2,1)==1 && size(point2,2)>=3 && size(point1,1)>0 && size(point1,2)>=3
    dist1 = sqrt((point1(:,1)-point2(1,1)).^2+(point1(:,2)-point2(1,2)).^2+(point1(:,3)-point2(1,3)).^2);
else
    dist1 = [];
end