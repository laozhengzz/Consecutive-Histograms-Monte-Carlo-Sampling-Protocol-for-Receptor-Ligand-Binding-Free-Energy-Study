function ring_atom = ring_judge(input_structure)

ring_atom = zeros(100,6);

ran=1;
for i = 1:size(input_structure,1)
    clear candidate_atom

    a = 1;
    candidate_atom(a,:) = input_structure(i,:);

    a = a+1;
    for j = 1:size(input_structure,1)
        if norm(input_structure(j,1:3)-input_structure(i,1:3))<=5 && i~=j
            candidate_atom(a,:) = input_structure(j,:);

            a = a+1;
        end
    end
    trig = 0;
    for k = 2:size(candidate_atom,1)
        if norm(candidate_atom(k,1:3) - candidate_atom(1,1:3))<1.8
            p_1 = candidate_atom(k,:);
            %trig = 1;
            for l = 2:size(candidate_atom,1)
                if norm(candidate_atom(l,1:3) - p_1(1,1:3))<1.8 && k ~=l
                    p_2 = candidate_atom(l,:);
                    if norm(candidate_atom(l,1:3) - candidate_atom(1,1:3))<1.8
                        ring_atom(ran,1) = input_structure(i,7);
                        ring_atom(ran,2) = candidate_atom(k,7);
                        ring_atom(ran,3) = candidate_atom(l,7);
                        ran=ran+1;
                        continue
                    end

                    for m = 2:size(candidate_atom,1)
                        if norm(candidate_atom(m,1:3) - p_2(1,1:3))<1.8 && l ~=m && k ~=m
                            p_3 = candidate_atom(m,:);
                            if norm(candidate_atom(m,1:3) - candidate_atom(1,1:3))<1.8
                                ring_atom(ran,1) = input_structure(i,7);
                                ring_atom(ran,2) = candidate_atom(k,7);
                                ring_atom(ran,3) = candidate_atom(l,7);
                                ring_atom(ran,4) = candidate_atom(m,7);
                                ran=ran+1;
                                continue
                            end



                            for n = 2:size(candidate_atom,1)
                                if norm(candidate_atom(n,1:3) - p_3(1,1:3))<1.8 && m ~=n && l ~=n && k ~=n
                                    p_4 = candidate_atom(n,:);
                                    if norm(candidate_atom(n,1:3) - candidate_atom(1,1:3))<1.8
                                        ring_atom(ran,1) = input_structure(i,7);
                                        ring_atom(ran,2) = candidate_atom(k,7);
                                        ring_atom(ran,3) = candidate_atom(l,7);
                                        ring_atom(ran,4) = candidate_atom(m,7);
                                        ring_atom(ran,5) = candidate_atom(n,7);
                                        ran=ran+1;
                                        continue
                                    end

                                    for p = 2:size(candidate_atom,1)
                                        if norm(candidate_atom(p,1:3) - p_4(1,1:3))<1.8 && n ~=p && m ~=p && l ~=p && k ~=p
                                            p_4 = candidate_atom(p,:);
                                            if norm(candidate_atom(p,1:3) - candidate_atom(1,1:3))<1.8
                                                ring_atom(ran,1) = input_structure(i,7);
                                                ring_atom(ran,2) = candidate_atom(k,7);
                                                ring_atom(ran,3) = candidate_atom(l,7);
                                                ring_atom(ran,4) = candidate_atom(m,7);
                                                ring_atom(ran,5) = candidate_atom(n,7);
                                                ring_atom(ran,6) = candidate_atom(p,7);
                                                ran=ran+1;
                                                continue
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


ring_atom(ring_atom(:,1) == 0,:)=[];

ring_atom = ring_atom';

for i = 1:size(ring_atom,2)
    ring_atom(:,i) = sort(ring_atom(:,i));
end
ring_atom = ring_atom';
ring_atom = unique(ring_atom, 'rows');



