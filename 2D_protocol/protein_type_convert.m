function protein_sol = protein_type_convert(protein_atom,AtomName,resName,protein_resseq)


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

protein_ring=zeros(1,length(protein_atom));
protein_ring_size=zeros(1,length(protein_atom));
for i=1:size(protein_atom,1)
    if (strncmp(AtomName{i,1}, hydrogen, 1)==1)|(strncmp(AtomName{i,1}, '1', 1)==1)|(strncmp(AtomName{i,1}, '2', 1)==1)|(strncmp(AtomName{i,1}, '3', 1)==1)|(strcmp(AtomName{i,1}, 'H')==1)
    %if (strcmp(AtomName{i,1}, 'H')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;    
    elseif (strcmp(AtomName{i,1}, 'HN1')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;     
    elseif (strcmp(AtomName{i,1}, 'HN2')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;    
    elseif (strcmp(AtomName{i,1}, 'HN3')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;        
        
    

    elseif (strcmp(AtomName{i,1}, 'O')==1) && (strncmp(resName{i,1}, 'HOH', 3)==1)% && (length(AtomName{i,1})==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.6546;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'NE', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=1;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NH1', 3)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N_pl3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NH2', 3)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N_pl3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=C_3Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    

    elseif (strncmp(AtomName{i,1}, 'OD1', 3)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=O_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'ND2', 3)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=N_pl3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
        
    elseif (strncmp(AtomName{i,1}, 'OD1', 3)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=O_co2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OD2', 3)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=O_co2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    %elseif (strncmp(AtomName{i,1}(1), 'S', 1)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
    elseif (strncmp(AtomName{i,1}, 'S', 1)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)    
        protein_name(i,1)=S_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(19,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.087;
        protein_backbone(i)=0;
        
        
        
    elseif (strncmp(AtomName{i,1}, 'OE1', 3)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=O_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=N_pl3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'OE1', 3)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=O_co2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OE2', 3)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=O_co2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=N_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(9,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
        
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=1;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=N_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(13,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
        
        
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=1;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=N_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(13,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
        
        
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=1;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
        
        
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CG1', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'NZ', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=N_4;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(11,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=1;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=C_3Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'S', 1)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=S_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(19,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    
    elseif (strcmp(AtomName{i,1}, 'CA')==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=C_3Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=C_3     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=C_3Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'N', 1)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=N_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(13,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=C_3Na ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OG', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=O_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=C_3Na ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OG1', 3)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=O_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'NE1', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_arNa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
    elseif (strncmp(AtomName{i,1}, 'CE3', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CZ2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CZ3', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CH2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=C_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=5;
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C_ar;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C_arNa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        protein_ring(i)=1;
        protein_ring_size(i)=6;
    elseif (strncmp(AtomName{i,1}, 'OH', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=O_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CG1', 3)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=C_3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'LPD', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'LPG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
        
    
    elseif (strncmp(AtomName{i,1}, 'OXT', 3)==1)
        protein_name(i,1)=O_2 ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
        
    elseif (strcmp(AtomName{i,1}, 'CA')==1)
        protein_name(i,1)=C_3Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1)
        protein_name(i,1)=C_3     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, 'C')==1)%&(length(AtomName{i,1})==1)
        protein_name(i,1)=C_2Na;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1)%&(length(AtomName{i,1})==1)
        protein_name(i,1)=O_2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1)%&(length(AtomName{i,1})==1)
        protein_name(i,1)=N_am;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0;
        protein_backbone(i)=1;
    elseif (strncmp(AtomName{i,1}, 'ZN', 2)==1) && (strncmp(resName{i,1}, 'ZN', 2)==1)
        protein_name(i,1)=ZN;
        protein_name(i,2)=ZN;
        protein_atom_radius(i)=1.2;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=2;
        protein_backbone(i)=0;
    end
        
    
    %protein_resseq(i)=PDBData1.Atom(1,i).resSeq;
    
    
end
protein_atom_num=[1:1:size(protein_atom,1)]';


protein_sol=[protein_atom,protein_name,protein_atom_radius',protein_atom_num,transpose(protein_ring),transpose(protein_ring_size),transpose(protein_atom_hb)];

protein_sol(protein_sol(:,4) ==0 | protein_sol(:,5) ==0,:)=[];
protein_sol(protein_sol(:,4)==35|(protein_sol(:,4)>=22&protein_sol(:,4)<=25),:)=[];

protein_sol(:,10) = 1:size(protein_sol,1);  
