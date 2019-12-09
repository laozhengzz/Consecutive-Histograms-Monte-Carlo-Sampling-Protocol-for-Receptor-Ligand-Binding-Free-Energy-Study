function an = ligandsig(ligand)
priority_atomtype_list=[8	10
    7	10
    9	10
    13	10
    12	10
    5	10
    11	10
    10	10
    16	5
    17	5
    18	5
    6	2
    32	2
    33	2
    34	2
    29	2
    30	2
    31	2
    26	2
    27	2
    28	2
    20	1
    15	1
    4	2
    3	1
    2	2
    1	1
    19	0
    14	0
    21	0
    ];

an=0;
for i = 1:size(ligand,1)
    for j =1:size(priority_atomtype_list( priority_atomtype_list(:,2)==10,:),1)
        if ligand(i,4) == priority_atomtype_list(j,1)
            %ligand_sig(an,1)=ligand(i,7);
            
            an=an+1;
            break
        end
    end
end

