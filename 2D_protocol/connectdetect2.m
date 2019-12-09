function connect_num = connectdetect2(atom1,atom2,R)

connect_num = 0;
%for i = 1:size(atom1,1)
    for j = 1:size(atom2,1)
        if norm(atom1(1,1:3)-atom2(j,1:3))<R && norm(atom1(1,1:3)-atom2(j,1:3))>1
            connect_num = connect_num+1;
        end
    end
%end