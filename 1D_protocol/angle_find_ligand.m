function angle_list = angle_find_ligand(input_structure)

%input_structure = protential_protein;
input_structure(:,10) = 0;
tl=1;
for i = 1:size(input_structure,1)
    clear candidate_atom

    a = 1;
    candidate_atom(a,:) = input_structure(i,:);
    a = a+1;
    for j = 1:size(input_structure,1)
        if sqrt((input_structure(j,1)-input_structure(i,1))^2 + (input_structure(j,2)-input_structure(i,2))^2 + (input_structure(j,3)-input_structure(i,3))^2 )<=4 && i~=j
            candidate_atom(a,:) = input_structure(j,:);
            a = a+1;
        end
    end
    trig = 0;
    for k = 2:size(candidate_atom,1)
        if (dist(candidate_atom(k,:),candidate_atom(1,:))<1.8 && (candidate_atom(k,4) ~= 20 && candidate_atom(k,4) ~= 21 && candidate_atom(1,4) ~= 20 && candidate_atom(1,4) ~= 21)) || (dist(candidate_atom(k,:),candidate_atom(1,:))<2.1 && (candidate_atom(k,4) == 20 || candidate_atom(k,4) == 21 || candidate_atom(1,4) == 20 || candidate_atom(1,4) == 21))
        
            p_1 = candidate_atom(k,:);
            trig = 0;
            for l = 2:size(candidate_atom,1)
                if ((dist(candidate_atom(l,:),p_1(1,:))<1.8 && (candidate_atom(l,4) ~= 20 && candidate_atom(l,4) ~= 21 && p_1(1,4) ~= 20 && p_1(1,4) ~= 21)) || (dist(candidate_atom(l,:),p_1(1,:))<2.1 && (candidate_atom(l,4) == 20 || candidate_atom(l,4) == 21 || p_1(1,4) == 20 || p_1(1,4) == 21))) && k ~=l
                    p_2 = candidate_atom(l,:);
                    angle_list(tl,1)=candidate_atom(1,4);
                    angle_list(tl,2)=p_1(1,4);
                    angle_list(tl,3)=p_2(1,4);
                    
                    angle_list(tl,4)=norm(candidate_atom(1,1:3) - p_2(1,1:3));
                    
                    angle_list(tl,5)=candidate_atom(1,7);
                    angle_list(tl,6)=p_2(1,7);
                    tl = tl+1;
                    trig = trig+1;
                    if trig ==3
                        break
                    end
                    

                end
            end
            %if trig ==1
            %    break
            %end

        end
    end


end

for n = 1:size(angle_list,1)
    if angle_list(n,1) > angle_list(n,3)
        temp_list = angle_list(n,:);
        angle_list(n,1) = temp_list(1,3);
        angle_list(n,2) = temp_list(1,2);
        angle_list(n,3) = temp_list(1,1);
        
        
        angle_list(n,4) = temp_list(1,4);
        angle_list(n,5) = temp_list(1,6);
        angle_list(n,6) = temp_list(1,5);
        
        
    end
end

angle_list = unique(angle_list,'rows');

















