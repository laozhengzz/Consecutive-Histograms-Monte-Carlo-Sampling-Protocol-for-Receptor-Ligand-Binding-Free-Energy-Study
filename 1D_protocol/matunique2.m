function data_table = matunique2(input_table,Col_num)

input_table = sortrows(input_table,[Col_num(1,1) Col_num(1,2)]);

mi = 1;

while mi < size(input_table,1)
    ft=0;
    for me = 1:size(input_table,1)
        if mi+me < size(input_table,1)
            
            if sum(abs(input_table(mi,:) - input_table(mi+me,:))) >= 10^-3
                input_table(mi+1:mi+me-1,:)=0;
                mi = mi+me;
                ft=0;
                break
            elseif sum(abs(input_table(mi,:) - input_table(mi+me,:))) < 10^-3
                ft=1;
                continue
            end
            
        elseif mi+me == size(input_table,1) && ft==1 && sum(abs(input_table(mi,:) - input_table(mi+me,:))) < 10^-3
            input_table(mi+1:mi+me,:)=0;
            mi = mi+me;
            ft=0;
            break
        elseif mi+me == size(input_table,1) && ft==1 && sum(abs(input_table(mi,:) - input_table(mi+me,:))) >= 10^-3
            input_table(mi+1:mi+me-1,:)=0;
            mi = mi+me;
            ft=0;
            break
            
        elseif mi+me == size(input_table,1) && ft==0
            mi = mi+me;
            ft=0;
            break
        elseif mi+me > size(input_table,1)
            mi = mi+me;
            break
            
        end
    end
    
end


input_table(input_table(:,1)==0,:)=[];

data_table = input_table;


