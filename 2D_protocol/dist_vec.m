function dist1 = dist_vec(point1,point2)
dist1 = [];
for i = 1:size(point1,1)
    dist1(i,1) = sqrt((point1(i,1)-point2(i,1)).^2+(point1(i,2)-point2(i,2)).^2+(point1(i,3)-point2(i,3)).^2);
end


