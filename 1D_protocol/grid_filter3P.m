function bin_s = grid_filter3P(grid_group1)
bin_s=[];
n_step=3;
%if size(grid_group1,1)>15
xmin_group1 = min(grid_group1(:,1))-0.001;
xmax_group1 = max(grid_group1(:,1))+0.001;

ymin_group1 = min(grid_group1(:,2))-0.001;
ymax_group1 = max(grid_group1(:,2))+0.001;

zmin_group1 = min(grid_group1(:,3))-0.001;
zmax_group1 = max(grid_group1(:,3))+0.001;

step_x=abs(xmin_group1-xmax_group1)/n_step;
step_y=abs(ymin_group1-ymax_group1)/n_step;
step_z=abs(zmin_group1-zmax_group1)/n_step;
sig_x = (xmin_group1:step_x:xmax_group1)';
sig_y = (ymin_group1:step_y:ymax_group1)';
sig_z = (zmin_group1:step_z:zmax_group1)';

bins=1;

for binx = 1:size(sig_x,1)-1
    for biny = 1:size(sig_y,1)-1
        for binz = 1:size(sig_z,1)-1
            bin_s1 = grid_group1((grid_group1(:,1)>sig_x(binx,1) & grid_group1(:,1)<=sig_x(binx+1,1)) & (grid_group1(:,2)>sig_y(biny,1) & grid_group1(:,2)<=sig_y(biny+1,1)) & (grid_group1(:,3)>sig_z(binz,1) & grid_group1(:,3)<=sig_z(binz+1,1)) ,:);
            if ~isempty(bin_s1)
                
                bin_s1 = sortrows(bin_s1,[4]);
                bin_s2 = bin_s1(1,:);
                %if size(bin_s2,1) ==1
                %    bin_s2 = bin_s2(1,:);
                %elseif size(bin_s2,1) >=2
                %    bin_s2 = bin_s2(1:2,:);
                %elseif size(bin_s2,1) >2
                %    bin_s2 = bin_s2(1:3,:);    
                %end
                
                bin_s(size(bin_s,1)+1:size(bin_s,1)+size(bin_s2,1),:)=bin_s2;
            end
        end
    end
end
%end
%if size(bin_s,1)>4
%    bin_s=bin_s(1:4,:);
%end