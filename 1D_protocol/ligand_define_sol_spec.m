function [ligand,name_str]= ligand_define_sol_spec(ligand_dir)

%run KECSA_GARFF_RMSD_zeroendAA
%run Solvation_Mat
%run BFE_Mat
%run PP_inter_FF


op_radi=[0
    3.325
3.73
3.65
3.635
3.66
3.155
2.895
2.85
2.955
2.785
2.84
2.845
2.735
2.82
2.745
3.015
2.76
3.78
];
%list=flipud(list);

C=1  ;
Cs=2 ;
CA=3 ;
CB=4 ;
CC=5 ;
CN=6 ;
CR=7 ;
CT=8 ;
CV=9 ;
CW=10;
H=11 ;
HO=12;
N=13;
N2=14;
N3=15;
NA=16;
NB=17;
O=18 ;
O2=19;
OH=20;
S=21;
SH=22;



C_3=1;
C_2=2;
C_1=3;
C_ar=4;
O_3=5;
O_3p=6;
O_2=7;
O_co2=8;
O_2v=9;
N_2=10;
N_am=11;
N_pl3=12;
N_4=13;
P=14;
F=15;
Cl=16;
Br=17;
I=18;
C_cat=19;
S_3=20;
S_o=21;
HNa=22;
HOa=23;
HSa=24;
H3 =25;
C_3Oa=26;
C_3Na=27;
C_3La=28;
C_2Oa=29;
C_2Na=30;
C_2La=31;
C_arOa=32;
C_arNa=33;
C_arLa=34;
HCa=35;

N_3=12;
N_1=10;
N_ar=10;
S_2=21;
S_o2=22;





K    = 29 ;    
MG   = 30  ;  
CAA   = 31	;
AL   = 32	;
MN   = 33  	;
FE   = 34    ;
CO   = 35    ;
NI   = 36    ;
CU   = 37    ;
ZN   = 38    ;
        
        
N_A=0;
D=1;
D2=2;
A=3;
DA=4;
        
ALA=101;
ARG=102;
ASN=103;
ASP=104;
CYS=105;
GLN=106;
GLU=107;
GLY=108;
HIS=109;
ILE=110;
LEU=111;
LYS=112;
MET=113;
PHE=114;
PRO=115;
SER=116;
THR=117;
TRP=118;
TYR=119;
VAL=120;



carbon='C';
oxygen='O';
nitrogen='N';
sulfer='S';
fluorin='F';
chlorine='Cl';
bromine='Br';
iodine='I';
phosphor='P';
hydrogen='H';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fp=fopen(strcat(home_dir,'GARF_MT_dG.txt'),'wt');%'a'??????a.txt,??????????


clear atom_conj_num atom name name_str ligand_atom_num atom_contact contact_type contact_name atom_contact_num mass;
fidin=fopen(ligand_dir(1,:));
tline=fgetl(fidin);
k=1;

while ~feof(fidin)
    if length(tline)>=60
        m=1;
        ct=0;
        while m<length(tline)
            if (strcmp(tline(m), ' ')~=1) && ct==0
                ct=1;
                m=m+1;
            elseif (strcmp(tline(m), ' ')==1) && ct==1
                ct=2;
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==2
                ct=3;
                m=m+1;
            elseif (strcmp(tline(m), ' ')==1) && ct==3
                ct=4;
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==4
                ct=5;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=6;
                        break
                    end
                end
                atom(k,1)=str2num(tline(m:m+mi-1));
                m=m+mi;
            elseif (strcmp(tline(m), ' ')==1) && ct==6
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==6
                
                ct=7;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=8;
                        break
                    end
                end
                atom(k,2)=str2num(tline(m:m+mi-1));
                m=m+mi;
            elseif (strcmp(tline(m), ' ')==1) && ct==8
                m=m+1;
                
            elseif (strcmp(tline(m), ' ')~=1) && ct==8
                
                ct=9;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=10;
                        break
                    end
                end
                atom(k,3)=str2num(tline(m:m+mi-1));
                m=m+mi;    
            elseif (strcmp(tline(m), ' ')==1) && ct==10
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==10
                name_str(k,1:5)=tline(m:m+4);
                
                ct=11;
                m=m+1;
            elseif ct==11
                break
            else
                m=m+1;
            end
        end
                
        k=k+1;
        tline=fgetl(fidin);
                
    elseif (length(tline)==13) && (strcmp(tline(1:13), '@<TRIPOS>BOND')==1)
        break
    elseif (length(tline)==21) && (strcmp(tline(1:21), '@<TRIPOS>SUBSTRUCTURE')==1)
        break
    else
        tline=fgetl(fidin);
    end
end

size_atom=size(atom);

size_name_str=size(name_str);

fclose(fidin);





%ligand initialization

clear name;  

size_name_str=size(name_str);

for i=1:1:size_name_str(1,1)
    o_pt=0;
    o2_pt=0;
    h_pt=0;
    c3_pt=0;
    c2_pt=0;
    car_pt=0;
    n_hn = 0;
    if strcmp(name_str(i,1:5), 'H    ')==1
        for j=1:size_atom(1,1)
            
            if strcmp(name_str(j,1), 'N')==1&norm(atom(j,1:3)-atom(i,1:3))<1.5
                name(i,1)=HNa;
                h_pt=1;
                break
                
            elseif strcmp(name_str(j,1), 'O')==1&norm(atom(j,1:3)-atom(i,1:3))<1.5
                name(i,1)=HOa;
                h_pt=1;
                break
                
            elseif strcmp(name_str(j,1), 'S')==1&norm(atom(j,1:3)-atom(i,1:3))<1.5
                name(i,1)=HSa;
                h_pt=1;
                break
                
            elseif (strcmp(name_str(j,1:5), 'N.4  ')==1|strcmp(name_str(j,1:5), 'C.cat')==1|strcmp(name_str(j,1:5), 'N_4  ')==1|strcmp(name_str(j,1:5), 'C_cat')==1)&norm(atom(j,1:3)-atom(i,1:3))<1.5
                name(i,1)=H3;
                h_pt=1;
                break
            end
        end
        
        if h_pt==0
            name(i,1)=H3;
        end
    elseif strcmp(name_str(i,1:5), 'HNa  ')==1
        name(i,1)=HNa;
    elseif strcmp(name_str(i,1:5), 'HOa  ')==1
        name(i,1)=HOa;
    elseif strcmp(name_str(i,1:5), 'HSa  ')==1
        name(i,1)=HSa;
    elseif strcmp(name_str(i,1:5), 'H3   ')==1
        name(i,1)=H3;
    elseif strcmp(name_str(i,1:5), 'C.3  ')==1 || strcmp(name_str(i,1:5), 'C_3  ')==1
        for j=1:size_atom(1,1)
            
            if (strcmp(name_str(j,1), 'N')==1 | strcmp(name_str(j,1), 'O')==1|strcmp(name_str(j,1), 'F')==1|strcmp(name_str(j,1:2), 'Cl')==1 ) & norm(atom(j,1:3)-atom(i,1:3))<1.9
                name(i,1)=C_3Na;
                c3_pt=1;
                break
            elseif (strcmp(name_str(j,1:2), 'Br')==1|strcmp(name_str(j,1), 'I')==1) & norm(atom(j,1:3)-atom(i,1:3))<2
                name(i,1)=C_3Na;
                c3_pt=1;
                break
            end
        end
        if c3_pt==0
            name(i,1)=C_3;
        end
    elseif strcmp(name_str(i,1:5), 'C_3Na')==1
        name(i,1)=C_3Na;
    elseif strcmp(name_str(i,1:5), 'C_3Oa')==1
        name(i,1)=C_3Na;
    elseif strcmp(name_str(i,1:5), 'C_3La')==1
        name(i,1)=C_3Na;
    elseif strcmp(name_str(i,1:5), 'C.2  ')==1 || strcmp(name_str(i,1:5), 'C_2  ')==1
        for j=1:size_atom(1,1)
            
            if (strcmp(name_str(j,1), 'N')==1 | strcmp(name_str(j,1), 'O')==1|strcmp(name_str(j,1), 'F')==1|strcmp(name_str(j,1:2), 'Cl')==1 ) &norm(atom(j,1:3)-atom(i,1:3))<1.9
                name(i,1)=C_2Na;
                c2_pt=1;
                break
            elseif (strcmp(name_str(j,1:2), 'Br')==1|strcmp(name_str(j,1), 'I')==1) & norm(atom(j,1:3)-atom(i,1:3))<2
                name(i,1)=C_2Na;
                c2_pt=1;
                break
            end
        end
        if c2_pt==0
            name(i,1)=C_2;
        end
        
    elseif strcmp(name_str(i,1:5), 'C_2Na')==1
        name(i,1)=C_2Na;
    elseif strcmp(name_str(i,1:5), 'C_2Oa')==1
        name(i,1)=C_2Na;
    elseif strcmp(name_str(i,1:5), 'C_2La')==1
        name(i,1)=C_2Na;
        
    elseif strcmp(name_str(i,1:5), 'C.1  ')==1 || strcmp(name_str(i,1:5), 'C_1  ')==1
        name(i,1)=C_1;
    elseif strcmp(name_str(i,1:5), 'C.ar ')==1 || strcmp(name_str(i,1:5), 'C_ar ')==1
        for j=1:size_atom(1,1)
            
            if (strcmp(name_str(j,1), 'N')==1 | strcmp(name_str(j,1), 'O')==1|strcmp(name_str(j,1), 'F')==1|strcmp(name_str(j,1:2), 'Cl')==1 ) &norm(atom(j,1:3)-atom(i,1:3))<1.9
                name(i,1)=C_arNa;
                car_pt=1;
                break
            elseif (strcmp(name_str(j,1:2), 'Br')==1|strcmp(name_str(j,1), 'I')==1) & norm(atom(j,1:3)-atom(i,1:3))<2
                name(i,1)=C_arNa;
                car_pt=1;
                break
            end
        end
        if car_pt==0
            name(i,1)=C_ar;
        end
        
    elseif strcmp(name_str(i,1:5), 'C_arN')==1
        name(i,1)=C_arNa;
    elseif strcmp(name_str(i,1:5), 'C_arO')==1
        name(i,1)=C_arNa;
    elseif strcmp(name_str(i,1:5), 'C_arL')==1
        name(i,1)=C_arNa;
    elseif strcmp(name_str(i,1:5), 'O.3  ')==1 || strcmp(name_str(i,1:5), 'O_3  ')==1
        for j=1:size_atom(1,1)
            if strcmp(name_str(j,1), 'H')==1&norm(atom(j,1:3)-atom(i,1:3))<1.5
                name(i,1)=O_3;
                o_pt=1;
                break
                
            end
        end
        if o_pt==0
            name(i,1)=O_3p;
            
        end
    elseif strcmp(name_str(i,1:5), 'O_3p ')==1
        name(i,1)=O_3p;
    elseif strcmp(name_str(i,1:5), 'O.2  ')==1 || strcmp(name_str(i,1:5), 'O_2  ')==1
        for j=1:size_atom(1,1)
            if (strcmp(name_str(j,1), 'P')==1|strcmp(name_str(j,1), 'S')==1)&norm(atom(j,1:3)-atom(i,1:3))<1.9
                name(i,1)=O_2v;
                o2_pt=1;
                break
                
            end
        end
        if o2_pt==0
            name(i,1)=O_2;
            
        end
    elseif strcmp(name_str(i,1:5), 'O_2v ')==1
        name(i,1)=O_2v;
        
        
        
    elseif strcmp(name_str(i,1:5), 'O.co2')==1 || strcmp(name_str(i,1:5), 'O_co2')==1
        name(i,1)=O_co2;
        
    elseif strcmp(name_str(i,1), 'N')==1
        for j=1:size_atom(1,1)
            if strcmp(name_str(j,1), 'H')==1&norm(atom(j,1:3)-atom(i,1:3))<1.5
                %name(i,1)=O_3;
                n_hn = n_hn+1;
                %break
            end
        end
        if n_hn == 0
            name(i,1)=N_2;
        elseif n_hn == 1
            name(i,1)=N_am;
        elseif n_hn == 2
            name(i,1)=N_pl3;
        elseif n_hn == 3
            name(i,1)=N_4;
        end
        
    elseif strcmp(name_str(i,1:5), 'N_2  ')==1
        name(i,1)=N_2;
    elseif strcmp(name_str(i,1:5), 'N_am ')==1
        name(i,1)=N_am;
    elseif strcmp(name_str(i,1:5), 'N_pl3')==1
        name(i,1)=N_pl3;
    elseif strcmp(name_str(i,1:5), 'N_4  ')==1
        name(i,1)=N_4;
        
    elseif strcmp(name_str(i,1:5), 'P.3  ')==1 || strcmp(name_str(i,1:5), 'P_3  ')==1 %&atom_charge(i,1)>=0.15&atom_charge(i,1)<=0.26
        name(i,1)=P;
        
    elseif strcmp(name_str(i,1:5), 'F    ')==1
        name(i,1)=F;
        
    elseif strcmp(name_str(i,1:5), 'Cl   ')==1%&atom_charge(i,1)>=0.15&atom_charge(i,1)<=0.26
        name(i,1)=Cl;
    elseif strcmp(name_str(i,1:5), 'Br   ')==1%&atom_charge(i,1)>=0.15&atom_charge(i,1)<=0.26
        name(i,1)=Br;
    elseif strcmp(name_str(i,1:5), 'I    ')==1%&atom_charge(i,1)>=0.15&atom_charge(i,1)<=0.26
        name(i,1)=I;
        
    elseif strcmp(name_str(i,1:5), 'C.cat')==1 || strcmp(name_str(i,1:5), 'C_cat')==1
        name(i,1)=C_cat;
        
    elseif strcmp(name_str(i,1:5), 'S.3  ')==1 || strcmp(name_str(i,1:5), 'S_3  ')==1
        name(i,1)=S_3;
    elseif strcmp(name_str(i,1:5), 'S.2  ')==1 || strcmp(name_str(i,1:5), 'S_2  ')==1
        name(i,1)=S_3;
    elseif strcmp(name_str(i,1:5), 'S.o  ')==1 || strcmp(name_str(i,1:5), 'S_o  ')==1
        name(i,1)=S_o;
    elseif strcmp(name_str(i,1:5), 'S.o2 ')==1 || strcmp(name_str(i,1:5), 'S_o2 ')==1
        name(i,1)=S_o;
    elseif strcmp(name_str(i,1:5), 'H    ')==1
        name(i,1)=HCa;
    elseif strcmp(name_str(i,1:5), 'HCa  ')==1
        name(i,1)=HCa;
    else
        name(i,1)=HCa;
    end
end
     

  

N_A=0;
D=1;
D2=2;
A=3;
DA=4;

clear atom_radius atom_hb acceptor_hb_angle donor_hb_angle atom_mass atom_charge

size_atom=size(atom);
for i=1:1:size_atom(1,1)
    if (name(i)==C_3) || (name(i)==C_3Oa) || (name(i)==C_3Na) || (name(i)==C_3La) 
        atom_radius(i)=op_radi(5,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=12;
        atom_charge(i)=0;
        
    elseif (name(i)==HCa) || (name(i)==HOa) || (name(i)==HSa) || (name(i)==HNa) || (name(i)==H3)
        atom_radius(i)=1;
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=1;
        atom_charge(i)=0;    
    elseif (name(i)==C_2) || (name(i)==C_2Oa) || (name(i)==C_2Na) || (name(i)==C_2La) 
        atom_radius(i)=op_radi(4,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=12;
        atom_charge(i)=0;    
        
    elseif (name(i)==C_1)
        atom_radius(i)=op_radi(3,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=12;
        atom_charge(i)=0;        
    elseif (name(i)==C_ar) || (name(i)==C_arOa) || (name(i)==C_arNa) || (name(i)==C_arLa) 
        atom_radius(i)=op_radi(6,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=12;
        atom_charge(i)=0;        
    elseif (name(i)==C_cat)
        atom_radius(i)=op_radi(6,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=12;
        atom_charge(i)=0;              
        
    elseif (name(i)==N_3)
        atom_radius(i)=op_radi(10,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=14;
        atom_charge(i)=0;
        
    elseif (name(i)==N_2)
        atom_radius(i)=op_radi(9,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=14;
        atom_charge(i)=0;    
        
    elseif (name(i)==N_ar)
        atom_radius(i)=op_radi(13,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=14;
        atom_charge(i)=0;        
        
    elseif (name(i)==N_4)
        atom_radius(i)=op_radi(11,1);
        atom_hb(i)=D;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=14;
        atom_charge(i)=1;
    elseif (name(i)==N_am)
        atom_radius(i)=op_radi(12,1);
        atom_hb(i)=D;
        atom_mass(i)=14;
        atom_charge(i)=0;
        
    elseif (name(i)==N_1)
        atom_radius(i)=op_radi(9,1);
        atom_hb(i)=N_A;
        atom_mass(i)=14;
        atom_charge(i)=0;
        
    elseif (name(i)==N_pl3)
        atom_radius(i)=op_radi(14,1);
        atom_hb(i)=D;
        atom_mass(i)=14;
        atom_charge(i)=0;    
        
    elseif ((name(i)==O_3))
        atom_radius(i)=op_radi(16,1);
        atom_hb(i)=DA;
        acceptor_hb_angle(i)=0.6083*pi;
        donor_hb_angle(i)=pi;
        atom_mass(i)=16;
        atom_charge(i)=0;
        
    elseif ((name(i)==O_3p))
        atom_radius(i)=op_radi(17,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=16;
        atom_charge(i)=0;
        
    elseif (name(i)==O_2)
        atom_radius(i)=op_radi(15,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=16;
        atom_charge(i)=0;    
        
    elseif (name(i)==O_2v)
        atom_radius(i)=op_radi(15,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0.6083*pi;
        donor_hb_angle(i)=pi;
        atom_mass(i)=16;
        atom_charge(i)=0;       
        
    elseif (name(i)==O_co2)
        atom_radius(i)=op_radi(18,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=pi;
        atom_mass(i)=16;
        atom_charge(i)=-1;    
        
        
        
        
    elseif name(i)==S_3
        atom_radius(i)=op_radi(19,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0.75*pi;
        donor_hb_angle(i)=0;
        atom_mass(i)=32;
        atom_charge(i)=0;
        
        

    elseif name(i)==S_2
        atom_radius(i)=op_radi(19,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0.75*pi;
        donor_hb_angle(i)=0;
        atom_mass(i)=32;
        atom_charge(i)=0;
    elseif name(i)==S_o
        atom_radius(i)=op_radi(19,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0.75*pi;
        donor_hb_angle(i)=0;
        atom_mass(i)=32;
        atom_charge(i)=0;
    elseif name(i)==S_o2
        atom_radius(i)=op_radi(19,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=32;
        atom_charge(i)=0;
    elseif name(i)==P
        atom_radius(i)=op_radi(19,1);
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=31;
        atom_charge(i)=0;
    elseif name(i)==F
        atom_radius(i)=op_radi(8,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=19;
        atom_charge(i)=0;
    elseif name(i)==Cl
        atom_radius(i)=op_radi(7,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=35.5;
        atom_charge(i)=0;
    elseif name(i)==Br
        atom_radius(i)=op_radi(2,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=80;
        atom_charge(i)=0;
    elseif name(i)==I
        atom_radius(i)=op_radi(2,1);
        atom_hb(i)=A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=127;
        atom_charge(i)=0;
    elseif name(i)==H
        atom_radius(i)=1.0;
        atom_hb(i)=N_A;
        acceptor_hb_angle(i)=0;
        donor_hb_angle(i)=0;
        atom_mass(i)=1;
        atom_charge(i)=0;
    end
end

total_atom_mass=sum(atom_mass);

clear neib_H

for i=1:1:size_atom(1,1)
    neib_H(i)=0;
    for j=1:1:size_atom(1,1)
    if ((name(i)==N_3) || (name(i)==N_2) || (name(i)==N_am) || (name(i)==N_pl3) || (name(i)==N_4)) && (name(j)==H) && (sqrt((atom(i,1)-atom(j,1))^2+(atom(i,2)-atom(j,2))^2+(atom(i,3)-atom(j,3))^2)<1.6)
        neib_H(i)=neib_H(i)+1;
    elseif ((name(i)==S_3) || (name(i)==S_2)) && (name(j)==H) && (sqrt((atom(i,1)-atom(j,1))^2+(atom(i,2)-atom(j,2))^2+(atom(i,3)-atom(j,3))^2)<1.6)
        neib_H(i)=neib_H(i)+1;
    end
    end
end



for i=1:1:size_atom(1,1)
    if (neib_H(i)==1) && ((name(i)==N_3) || (name(i)==N_2) || (name(i)==N_am) || (name(i)==N_pl3) || (name(i)==N_4))
        atom_hb(i)=D;
        donor_hb_angle(i)=pi;
    elseif (neib_H(i)==2) && ((name(i)==N_3) || (name(i)==N_2) || (name(i)==N_am) || (name(i)==N_pl3) || (name(i)==N_4))
        atom_hb(i)=D2;
        donor_hb_angle(i)=pi;
    elseif (neib_H(i)==1) && ((name(i)==S_3) || (name(i)==S_2))
        atom_hb(i)=DA;
        donor_hb_angle(i)=pi;
        acceptor_hb_angle(i)=0.6083*pi;
    end
end

       

for i=1:1:size_atom(1,1)
    for j=1:1:size_atom(1,1)
    if (name(i)==O_3p) && (name(j)==H) && (sqrt((atom(i,1)-atom(j,1))^2+(atom(i,2)-atom(j,2))^2+(atom(i,3)-atom(j,3))^2)<1.5)
        atom_radius(i)=op_radi(16,1);
       atom_hb(i)=DA;
       donor_hb_angle(i)=pi;
       acceptor_hb_angle(i)=0.6083*pi;
       name(i)=O_3;
       break
    end
    end
end
        
for i=1:1:size_atom(1,1)
    for j=1:1:size_atom(1,1)
    if (name(i)==O_2) && (name(j)==P|name(j)==S_o|name(j)==S_o2) && (sqrt((atom(i,1)-atom(j,1))^2+(atom(i,2)-atom(j,2))^2+(atom(i,3)-atom(j,3))^2)<1.9)
        atom_radius(i)=op_radi(15,1);
       atom_hb(i)=A;
       donor_hb_angle(i)=pi;
       acceptor_hb_angle(i)=0.6083*pi;
       name(i)=O_2v;
       break
    end
    end
end
        
   
        

clear ligand ligand_atom_num ligand_H
ligand=[atom,name,transpose(atom_mass),transpose(atom_radius),transpose(atom_hb),transpose(atom_hb),transpose(atom_hb),transpose(atom_hb)];

R1=2.0;
R2=2.1;

%for jz = 1:size(ligand,1)
%    if ligand(jz,4) ~= 1 && ligand(jz,4) ~= 27
%        continue
%    elseif ligand(jz,4) == 1 || ligand(jz,4) == 27
%        
%        connect_num = connectdetect2(ligand(jz,:),ligand,R1);
%        if connect_num == 3 && ligand(jz,4) == 1
%            ligand(jz,4) = 4;
%        elseif connect_num == 3 && ligand(jz,4) == 27
%            
%            ligand(jz,4) = 33;
%        end
%    end
%end


ligand(ligand(:,4)>=22 & ligand(:,4)<=25,:)=[];

size_ligand=size(ligand);

ligand_atom_num=[1:1:size_ligand(1,1)];

ligand(:,7)=transpose(ligand_atom_num);
ring_atom = ring_judge(ligand);

C2ar_atom = ligand(ligand(:,4)==2|ligand(:,4)==4|ligand(:,4)==30|ligand(:,4)==33,:);
R1 = 1.8;
for i = 1:size(ligand,1)
    lm = 0;
    for j = 1:size(ring_atom,1)
        for k = 1:size(ring_atom,2)
            if ligand(i,7) == ring_atom(j,k) && ligand(i,4) == 2 && size(ring_atom(j,ring_atom(j,:)>0),2) == 6
                connect_num = connectdetect2(ligand(i,:),C2ar_atom,R1);
                if connect_num>=2
                    
                    ligand(i,4) = 4;
                    lm = 1;
                    break
                end
            end
            if lm == 1
                break
            end
        end
    end
end

N_atom = ligand(ligand(:,4)>=10 & ligand(:,4)<=13,:);
O_atom = ligand(ligand(:,4)>=5 & ligand(:,4)<=9,:);
O2_atom = ligand(ligand(:,4)>=7 & ligand(:,4)<=9,:);
Hlg_atom = ligand(ligand(:,4)>=15 & ligand(:,4)<=18,:);

Polar_small_atom = [N_atom
    O_atom];

R1=2.0;
R2=2.1;
C2X_atom = ligand(ligand(:,4) == 30,:);

for jz=1:1:size(ligand,1)
    ring_mark = 0;
    ring_size = 0;
    for r1 = 1:size(ring_atom,1)
        for r2 = 1:size(ring_atom,2)
            if ring_atom(r1,r2) == ligand(jz,7);
                ring_mark = 1;
                ring_row = ring_atom(r1,:);
                ring_size = size(ring_row(:,ring_row(1,:)~=0),2);
                break
            end
        end
        if ring_mark == 1
            break
        end
    end
    ligand(jz,8) = ring_mark;
    ligand(jz,9) = ring_size;
    
    
    connect_num1 = connectdetect2(ligand(jz,:),Polar_small_atom,R1);
    connect_num2 = connectdetect2(ligand(jz,:),Hlg_atom,R2);
    %connect_num1 = connectdetect2(ligand(jz,:),O2_atom,R1);
    connect_num3 = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    if ligand(jz,4)==30 && (connect_num1 + connect_num2 + connect_num3)>=2 && ring_mark == 1
        
        ligand(jz,4)=33;
    end
end









