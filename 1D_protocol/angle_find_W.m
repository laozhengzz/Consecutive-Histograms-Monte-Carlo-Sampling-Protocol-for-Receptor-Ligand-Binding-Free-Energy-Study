function angle_mark = angle_find_W(atom,mid_structure,input_structure,R1,R2)

angle_mark = 0;
input_structure(:,10) = 0;
tl=1;

candidate_atom=[];

a = 1;

for j = 1:size(input_structure,1)
    if sqrt((input_structure(j,1)-atom(1,1))^2 + (input_structure(j,2)-atom(1,2))^2 + (input_structure(j,3)-atom(1,3))^2 )<=R2 && sqrt((input_structure(j,1)-atom(1,1))^2 + (input_structure(j,2)-atom(1,2))^2 + (input_structure(j,3)-atom(1,3))^2 )>0
        candidate_atom(a,:) = input_structure(j,:);
        a = a+1;
    end
end


mid_candidate_atom=[];

a = 1;

for j = 1:size(mid_structure,1)
    if sqrt((mid_structure(j,1)-atom(1,1))^2 + (mid_structure(j,2)-atom(1,2))^2 + (mid_structure(j,3)-atom(1,3))^2 )<=R2 && sqrt((mid_structure(j,1)-atom(1,1))^2 + (mid_structure(j,2)-atom(1,2))^2 + (mid_structure(j,3)-atom(1,3))^2 )>0
        mid_candidate_atom(a,:) = mid_structure(j,:);
        a = a+1;
    end
end







trig = 0;
for k = 1:size(mid_candidate_atom,1)
    if sqrt((mid_candidate_atom(k,1)-atom(1,1))^2 + (mid_candidate_atom(k,2)-atom(1,2))^2 + (mid_candidate_atom(k,3)-atom(1,3))^2 )<=R1
        p_1 = mid_candidate_atom(k,:);
        trig = 0;
        for l = 1:size(candidate_atom,1)
            if sqrt((candidate_atom(l,1)-p_1(1,1))^2 + (candidate_atom(l,2)-p_1(1,2))^2 + (candidate_atom(l,3)-p_1(1,3))^2 )<=R1
                
                angle_mark = angle_mark+1;
                tl = tl+1;
                trig = trig+1;
                if trig ==3
                    break
                end
                
                
            end
        end
        
    end
end
















