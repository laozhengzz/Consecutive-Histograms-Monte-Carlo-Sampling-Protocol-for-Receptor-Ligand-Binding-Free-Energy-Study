function [protein,protein_list_num] = protein_define_final_wwaterclose(protein_dir)

%clear all
%protein_dir='E:\Affinity_Benchmark\1VFB_r_b.pdb';
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

%original_dir = 'E:\SiMTile159\'
%list_ori=dir(original_dir);
%size_list_ori=size(list_ori);

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


%protein atoms initialization
clear PDBData1 protein_atom1 protein_atom2 protein_atom3 protein_atomname protein_resseq protein_occupancy protein_atom_f protein_atom


fidin=fopen(protein_dir(1,:));
tline=fgetl(fidin);
k=1;
num=1;
clear AtomName resName resSeq chainID protein_atom
while ~feof(fidin)
    
        

    if length(tline)>6 && ((strcmp(tline(1:4), 'ATOM')==1)||(strcmp(tline(1:6), 'HETATM')==1))
        m=1;
        ct=0;
        while m<length(tline)
            if (strcmp(tline(m), ' ')==1) && ct==0
                ct=1;
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==1
                ct=2;
                m=m+1;
            elseif (strcmp(tline(m), ' ')==1) && ct==2
                ct=3;
                m=m+1;
           
            elseif (strcmp(tline(m), ' ')~=1) && ct==3
                ct=4;
                for mi=1:3
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=5;
                        AtomName{k,1}=tline(m:m+mi-1);
                        m=m+mi;
                        break
                    end
                end
                if ct ==4
                    AtomName{k,1}=tline(m:m+mi);
                    ct=5;
                    m=m+mi+1;
                end
                %cutA(1,2)=m+mi-1;
                
                %m=m+mi;
            elseif (strcmp(tline(m), ' ')==1) && ct==5
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==5
                ct=6;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=7;
                        break
                    end
                end
                RN_test = tline(m:m+mi-1);
                if size(RN_test,2)>3
                    RN_test=RN_test(end-2:end);
                    
                end
                
                resName{k,:}=RN_test;
                
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==7
                m=m+3;   
                ct=9;
            %%%%%%%%%%%%%%%%%chain ID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %elseif (strcmp(tline(m), ' ')~=1) && ct==7    
            %    ct=8;
            %    for mi=1:10
            %        if (strcmp(tline(m+mi), ' ')==1)
            %            ct=9;
            %            break
            %        end
            %    end
            %    chainID{k,1}=tline(m:m+mi-1);
            %    m=m+mi;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (strcmp(tline(m), ' ')==1) && ct==9
                m=m+1;    
            elseif (strcmp(tline(m), ' ')~=1) && ct==9    
                ct=10;    
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=11;
                        break
                    end
                end
                str1=tline(m:m+mi-1);
                resName1(k,1)=str2num(str1(regexp(str1,'\d')));
                
                %resSeq(k,1)=str2num(str1(regexp(str1,'\d')));
                m=m+mi;
                
                
            elseif (strcmp(tline(m), ' ')==1) && ct==11
                m=m+1;    
            elseif (strcmp(tline(m), ' ')~=1) && ct==11
                ct=12;    
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=13;
                        break
                    end
                end
                protein_atom(k,1)=str2num(tline(m:m+mi-1));
                m=m+mi;    
                
            elseif (strcmp(tline(m), ' ')==1) && ct==13
                m=m+1;    
            elseif (strcmp(tline(m), ' ')~=1) && ct==13
                ct=14;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=15;
                        break
                    end
                end
                protein_atom(k,2)=str2num(tline(m:m+mi-1));
                m=m+mi;    
                
            elseif (strcmp(tline(m), ' ')==1) && ct==15
                m=m+1;    
            elseif (strcmp(tline(m), ' ')~=1) && ct==15
                ct=16;
                for mi=1:10
                    if length(tline) == (m+mi);
                        protein_atom(k,3)=str2num(tline(m:m+mi));
                        ct=17;
                        break
                        
                    elseif (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        protein_atom(k,3)=str2num(tline(m:m+mi-1));
                        ct=17;
                        break
                    end
                end
                
                m=m+mi;
            elseif ct==17
                break
            else
                m=m+1;
            end
        end
        
        k=k+1;
        tline=fgetl(fidin);
        
    elseif (length(tline)==3) && (strcmp(tline(1:3), 'END')==1)
        break
    else
        tline=fgetl(fidin);
    end
end



if (length(tline)>6)
    if (strcmp(tline(1:3), 'END')~=1) && ((strcmp(tline(1:4), 'ATOM')==1)||(strcmp(tline(1:6), 'HETATM')==1))
        m=1;
        ct=0;
        while m<length(tline)
            if (strcmp(tline(m), ' ')==1) && ct==0
                ct=1;
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==1
                ct=2;
                m=m+1;
            elseif (strcmp(tline(m), ' ')==1) && ct==2
                ct=3;
                m=m+1;
                
            elseif (strcmp(tline(m), ' ')~=1) && ct==3
                ct=4;
                for mi=1:3
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=5;
                        AtomName{k,1}=tline(m:m+mi-1);
                        break
                    end
                end
                if ct ==4
                    AtomName{k,1}=tline(m:m+mi);
                    ct=5;
                end
                %cutA(1,2)=m+mi-1;
                
                m=m+mi;
            elseif (strcmp(tline(m), ' ')==1) && ct==5
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==5
                ct=6;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=7;
                        break
                    end
                end
                RN_test = tline(m:m+mi-1);
                if size(RN_test,2)>3
                    RN_test=RN_test(end-2:end);
                    
                end
                
                resName{k,:}=RN_test;
                
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==7
                m=m+3;
                ct=9;
                %%%%%%%%%%%%%%%%%chain ID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %elseif (strcmp(tline(m), ' ')~=1) && ct==7
                %    ct=8;
                %    for mi=1:10
                %        if (strcmp(tline(m+mi), ' ')==1)
                %            ct=9;
                %            break
                %        end
                %    end
                %    chainID{k,1}=tline(m:m+mi-1);
                %    m=m+mi;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (strcmp(tline(m), ' ')==1) && ct==9
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==9
                ct=10;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1)
                        ct=11;
                        break
                    end
                end
                str1=tline(m:m+mi-1);
                resName1(k,1)=str2num(str1);
                
                resSeq(k,1)=str2num(str1(regexp(str1,'\d')));
                m=m+mi;
                
                
            elseif (strcmp(tline(m), ' ')==1) && ct==11
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==11
                ct=12;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=13;
                        break
                    end
                end
                protein_atom(k,1)=str2num(tline(m:m+mi-1));
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==13
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==13
                ct=14;
                for mi=1:10
                    if (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        ct=15;
                        break
                    end
                end
                protein_atom(k,2)=str2num(tline(m:m+mi-1));
                m=m+mi;
                
            elseif (strcmp(tline(m), ' ')==1) && ct==15
                m=m+1;
            elseif (strcmp(tline(m), ' ')~=1) && ct==15
                ct=16;
                for mi=1:10
                    if length(tline) == (m+mi);
                        protein_atom(k,3)=str2num(tline(m:m+mi));
                        ct=17;
                        break
                        
                    elseif (strcmp(tline(m+mi), ' ')==1) || (strcmp(tline(m+mi), '-')==1)
                        protein_atom(k,3)=str2num(tline(m:m+mi-1));
                        ct=17;
                        break
                    end
                end
                %protein_atom(k,3)=str2num(tline(m:m+mi-1));
                m=m+mi;
            elseif ct==17
                break
            else
                m=m+1;
            end
        end
        
        %k=k+1;
        %tline=fgetl(fidin);
    end
end





fclose(fidin);
        
for i = 1:size(resName,1)
    if strcmp(resName{i,1}, 'ALA')==1
        resSeq(i,1) = 101;
    elseif strcmp(resName{i,1}, 'ARG')==1
        resSeq(i,1) = 102;
    elseif strcmp(resName{i,1}, 'ASN')==1 || strcmp(resName{i,1}, 'ASX')==1
        resSeq(i,1) = 103;
    elseif strcmp(resName{i,1}, 'ASP')==1
        resSeq(i,1) = 104;
    elseif strcmp(resName{i,1}, 'CYS')==1
        resSeq(i,1) = 105;
    elseif strcmp(resName{i,1}, 'GLN')==1
        resSeq(i,1) = 106;
    elseif strcmp(resName{i,1}, 'GLU')==1 || strcmp(resName{i,1}, 'GLX')==1
        resSeq(i,1) = 107;
    elseif strcmp(resName{i,1}, 'GLY')==1
        resSeq(i,1) = 108;
    elseif strcmp(resName{i,1}, 'HIS')==1
        resSeq(i,1) = 109;
    elseif strcmp(resName{i,1}, 'HID')==1
        resSeq(i,1) = 109.1;
    elseif strcmp(resName{i,1}, 'HIE')==1
        resSeq(i,1) = 109.2;
    elseif strcmp(resName{i,1}, 'HIP')==1
        resSeq(i,1) = 109.3;    
    elseif strcmp(resName{i,1}, 'ILE')==1
        resSeq(i,1) = 110;
    elseif strcmp(resName{i,1}, 'LEU')==1
        resSeq(i,1) = 111;
    elseif strcmp(resName{i,1}, 'LYS')==1
        resSeq(i,1) = 112;
    elseif strcmp(resName{i,1}, 'MET')==1
        resSeq(i,1) = 113;
    elseif strcmp(resName{i,1}, 'PHE')==1
        resSeq(i,1) = 114;
    elseif strcmp(resName{i,1}, 'PRO')==1
        resSeq(i,1) = 115;
    elseif strcmp(resName{i,1}, 'SER')==1
        resSeq(i,1) = 116;
    elseif strcmp(resName{i,1}, 'THR')==1
        resSeq(i,1) = 117;
    elseif strcmp(resName{i,1}, 'TRP')==1
        resSeq(i,1) = 118;
    elseif strcmp(resName{i,1}, 'TYR')==1
        resSeq(i,1) = 119;
    elseif strcmp(resName{i,1}, 'VAL')==1
        resSeq(i,1) = 120;
    else
        resSeq(i,1) = 121;
    end
end
        
        
        
        
        clear protein_name protein_atom_radius protein_atom_hb protein_acceptor_hb_angle protein_donor_hb_angle protein_atom_charge protein_backbone
        
        protein_resseq=resSeq;
        size_protein_atom=size(protein_atom);
for i=1:size_protein_atom(1,1)
    %if (strncmp(AtomName{i,1}(1), hydrogen, 1)==1)|(strncmp(AtomName{i,1}(1), '1', 1)==1)|(strncmp(AtomName{i,1}(1), '2', 1)==1)|(strncmp(AtomName{i,1}(1), '3', 1)==1)
    %    protein_name(i,1)=H;
    %    protein_name(i,2)=protein_resseq(i);
    %    protein_atom_radius(i)=1;
    %    protein_atom_hb(i)=N_A;
    %    protein_acceptor_hb_angle(i)=0;
    %    protein_donor_hb_angle(i)=0;
    %    protein_atom_charge(i)=0;
    %    protein_backbone(i)=1;
    if (strcmp(AtomName{i,1}, 'H')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;    
    elseif (strcmp(AtomName{i,1}, 'HN1')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;     
    elseif (strcmp(AtomName{i,1}, 'HN2')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;    
    elseif (strcmp(AtomName{i,1}, 'HN3')==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=H;
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
        
        
    
    elseif (strcmp(AtomName{i,1}, 'C')==1)% && (length(AtomName{i,1})==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.623;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'O')==1)% && (length(AtomName{i,1})==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.574;
        protein_backbone(i)=1;
    elseif (strcmp(AtomName{i,1}, 'N')==1)% && (length(AtomName{i,1})==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.452;
        protein_backbone(i)=1;

        
    
    elseif (strncmp(AtomName{i,1}, 'NE', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.752;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HE', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.83;
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=0;    
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2637;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0007;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.039;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0486;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.8076;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NH1', 3)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.8627;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NH2', 3)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=N2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.8627;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, '1HH2')==1 || strcmp(AtomName{i,1}, '2HH2')==1 || strncmp(AtomName{i,1}, '1HH1',4)==1 || strncmp(AtomName{i,1}, '2HH1',4)==1 || strncmp(AtomName{i,1}, 'HH1',3)==1 || strncmp(AtomName{i,1}, 'HH2',3)==1 || strncmp(AtomName{i,1}, 'HH',2)==1) && (strncmp(resName{i,1}, 'ARG', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.83;
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4196;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0143;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2041;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OD1', 3)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5931;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'ND2', 3)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.9191;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, '1HD2')==1 || strcmp(AtomName{i,1}, '2HD2')==1 || strncmp(AtomName{i,1}, '1HD1',4)==1 || strncmp(AtomName{i,1}, '2HD1',4)==1 || strncmp(AtomName{i,1}, 'HD1',3)==1 || strncmp(AtomName{i,1}, 'HD2',3)==1 || strncmp(AtomName{i,1}, 'HD',2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.83;
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4196;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASN', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7130;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0143;
        protein_backbone(i)=1;     
        resName{i,1} = 'ASN';
        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2041;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
    elseif (strncmp(AtomName{i,1}, 'OD1', 3)==1 || strncmp(AtomName{i,1}, 'XD1', 3)==1 ) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5931;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
        AtomName{i,1} = 'OD1';
    elseif (strncmp(AtomName{i,1}, 'ND2', 3)==1 || strncmp(AtomName{i,1}, 'XD2', 3)==1 ) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.9191;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
        AtomName{i,1} = 'ND2';
    elseif (strcmp(AtomName{i,1}, '1HD2')==1 || strcmp(AtomName{i,1}, '2HD2')==1 || strncmp(AtomName{i,1}, '1HD1',4)==1 || strncmp(AtomName{i,1}, '2HD1',4)==1 || strncmp(AtomName{i,1}, 'HD1',3)==1 || strncmp(AtomName{i,1}, 'HD2',3)==1 || strncmp(AtomName{i,1}, 'HD',2)==1) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.83;
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4196;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASX', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7130;
        protein_backbone(i)=0;
        resName{i,1} = 'ASN';
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0381;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0303;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OD1', 3)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8014;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OD2', 3)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8014;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'ASP', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.7994;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0213;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1231;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}(1), 'S', 1)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=SH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(19,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.312;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}(1), 'HS', 2)==1) && (strncmp(resName{i,1}, 'CYS', 3)==1)
        protein_name(i,1)=HS;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=2.09;
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.193;
        protein_backbone(i)=0;    
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0031;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0036;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OE1', 3)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=O;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.6086;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=N;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(12,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.9407;
        protein_backbone(i)=0;
    elseif (strcmp(AtomName{i,1}, '1HE2')==1 || strcmp(AtomName{i,1}, '2HE2')==1 || strncmp(AtomName{i,1}, '1HE1',4)==1 || strncmp(AtomName{i,1}, '2HE1',4)==1 || strncmp(AtomName{i,1}, 'HE1',3)==1 || strncmp(AtomName{i,1}, 'HE2',3)==1 || strncmp(AtomName{i,1}, 'HE',2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.83;
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4251;
        protein_backbone(i)=0;    
        
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0645;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLN', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.6951;
        protein_backbone(i)=0;
        
        
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0397;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.056;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OE1', 3)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8188;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OE2', 3)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8188;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0136;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLU', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.8054;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0397;
        protein_backbone(i)=1;
        resName{i,1} = 'GLU';
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.056;
        protein_backbone(i)=0;
        resName{i,1} = 'GLU';
    elseif (strncmp(AtomName{i,1}, 'OE1', 3)==1 || strncmp(AtomName{i,1}, 'XE1', 3)==1) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8188;
        protein_backbone(i)=0;
        resName{i,1} = 'GLU';
        AtomName{i,1} = 'OE1';
    elseif (strncmp(AtomName{i,1}, 'OE2', 3)==1 || strncmp(AtomName{i,1}, 'XE2', 3)==1) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=O2;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(18,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.8188;
        protein_backbone(i)=0;
        resName{i,1} = 'GLU';
        AtomName{i,1} = 'OE2';
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0136;
        protein_backbone(i)=0;
        resName{i,1} = 'GLU';
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'GLX', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.8054;
        protein_backbone(i)=0;
        resName{i,1} = 'GLU';
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0581;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0074;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=NB;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(9,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.5432;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=CW;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2207;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=CR;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1635;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2795;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HIE', 3)==1)
        protein_name(i,1)=CC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1868;
        protein_backbone(i)=0;
   
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0188;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0462;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3811;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=CV;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1292;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=CR;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2057;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=NB;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(9,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5727;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HID', 3)==1)
        protein_name(i,1)=CC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0266;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0188;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0462;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3811;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=CV;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1292;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=CR;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2057;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=NB;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(9,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.5727;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HIS', 3)==1)
        protein_name(i,1)=CC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0266;
        protein_backbone(i)=0;
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1354;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0414;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'ND1', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.1513;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=CV;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1141;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=CR;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.017;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE2', 3)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1718;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'HIP', 3)==1)
        protein_name(i,1)=CC;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0012;
        protein_backbone(i)=0;
        
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0597;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1303;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG1', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.043;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.3204;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'ILE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.066;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0518;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1103;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.4121;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.4121;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'LEU', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3531;
        protein_backbone(i)=0;
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.24;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0094;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NZ', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=N3;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(11,1);
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3854;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HZ', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=H;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.83;
        protein_atom_hb(i)=D2;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.34;
        protein_backbone(i)=0;    
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0187;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.05;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1) && (strncmp(resName{i,1}, 'LYS', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0143;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.024;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0342;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0018;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}(1), 'S', 1)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=S;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(19,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2737;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'MET', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0536;
        protein_backbone(i)=0;
    
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0024;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0342;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0118;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1256;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1256;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1704;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1704;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'PHE', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1072;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.026;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.007;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0192;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'PRO', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0189;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0249;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2117;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OG', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.6546;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HG', 2)==1) && (strncmp(resName{i,1}, 'SER', 3)==1)
        protein_name(i,1)=HO;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4275;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0389;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3654;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OG1', 3)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.6761;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HG', 2)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=HO;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.4102;
        protein_backbone(i)=0;    
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'THR', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2438;
        protein_backbone(i)=0;
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0275;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.005;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CW;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1638;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CB;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.1243;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'NE1', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=NA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(14,1);
        protein_atom_hb(i)=D;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.3418;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CN;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.138;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE3', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2387;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2601;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ3', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1972;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CH2', 3)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1134;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'TRP', 3)==1)
        protein_name(i,1)=Cs;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1415;
        protein_backbone(i)=0;
        
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0014;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0152;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0011;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD1', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1906;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CD2', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.1906;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE1', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2341;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CE2', 3)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=CA;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(6,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.2341;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CZ', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=C;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(4,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.3226;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OH', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=OH;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(16,1);
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=-0.5579;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'HH', 2)==1) && (strncmp(resName{i,1}, 'TYR', 3)==1)
        protein_name(i,1)=HO;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1.74;
        protein_atom_hb(i)=DA;
        protein_acceptor_hb_angle(i)=0.6083*pi;
        protein_donor_hb_angle(i)=pi;
        protein_atom_charge(i)=0.3992;
        protein_backbone(i)=0;    
        
        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0875;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.2985;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG1', 3)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.3192;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'CG2', 3)==1) && (strncmp(resName{i,1}, 'VAL', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.3192;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'GLY', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.0252;
        protein_backbone(i)=1;        
        
    elseif (strncmp(AtomName{i,1}, 'CA', 2)==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0603;
        protein_backbone(i)=1;        
    elseif (strncmp(AtomName{i,1}, 'CB', 2)==1) && (strncmp(resName{i,1}, 'ALA', 3)==1)
        protein_name(i,1)=CT     ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.0603;
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
        
    

    elseif (strncmp(AtomName{i,1}, 'CE', 2)==1)
        protein_name(i,1)=CT;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(5,1);
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.045;
        protein_backbone(i)=0;
    elseif (strncmp(AtomName{i,1}, 'OXT', 3)==1)
        protein_name(i,1)=O ;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=op_radi(15,1);
        protein_atom_hb(i)=A;
        protein_acceptor_hb_angle(i)=0.75*pi;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=-0.574;
        protein_backbone(i)=1;
        
        
    elseif (strncmp(AtomName{i,1}, 'ZN', 2)==1)
        protein_name(i,1)=ZN;
        protein_name(i,2)=ZN;
        protein_atom_radius(i)=1.2;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=2;
        protein_backbone(i)=0;
        
    elseif (strncmp(AtomName{i,1}, 'H',1)==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;     
        
    elseif (strncmp(AtomName{i,1}, '1H',2)==1 || strncmp(AtomName{i,1}, '2H',2)==1 || strncmp(AtomName{i,1}, '3H',2)==1)% && (strncmp(PDBData1.Atom(1,i-1).AtomName, 'N',1)==1)
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;      
        
    else
        protein_name(i,1)=HCa;
        protein_name(i,2)=protein_resseq(i);
        protein_atom_radius(i)=1;
        protein_atom_hb(i)=N_A;
        protein_acceptor_hb_angle(i)=0;
        protein_donor_hb_angle(i)=0;
        protein_atom_charge(i)=0.324;
        protein_backbone(i)=1;      
    end
        
    %protein_resseq(i)=PDBData1.Atom(1,i).resSeq;
    
end

clear protential_protein protein_atom_num protein_H;
protential_protein=[protein_atom,protein_name,protein_atom_radius',transpose(protein_atom_hb),transpose(protein_backbone),resName1,transpose(protein_atom_hb)];

%for i = 1:size(protential_protein,1)
%    if protential_protein(i,4) ==0 || protential_protein(i,5) ==0 || protential_protein(i,4)==35 || protential_protein(i,4)==11 || protential_protein(i,4)==12
    %if protential_protein(i,4)==35 || protential_protein(i,4)==11 || protential_protein(i,4)==12
        AtomName((protential_protein(:,4) ==0 | protential_protein(:,5) ==0 | protential_protein(:,4)==35 | protential_protein(:,4)==11 | protential_protein(:,4)==12),:)=[];
        resName((protential_protein(:,4) ==0 | protential_protein(:,5) ==0 | protential_protein(:,4)==35 | protential_protein(:,4)==11 | protential_protein(:,4)==12),:)=[];
%    end
%end


protein_Namelist = [AtomName resName];

protein_list_num = protein_list_final_AF_format(protein_Namelist);


protential_protein(protential_protein(:,4) ==0 | protential_protein(:,5) ==0,:)=[];

size_protential_protein=size(protential_protein);

protein_atom_num=[1:1:size_protential_protein(1,1)];

protential_protein(:,7)=transpose(protein_atom_num);
protein_H=protential_protein(protential_protein(:,4)==35|protential_protein(:,4)==11|protential_protein(:,4)==12,:);

size_protein_H=size(protein_H);

protential_protein(protential_protein(:,4)==35|protential_protein(:,4)==11|protential_protein(:,4)==12,:)=[];

protein=protential_protein;


protein(:,10) = 1:size(protein,1);  


protein_water_all = protein(protein(:,4) == 20 & protein(:,5)==121,:);

protein(protein(:,4) == 20 & protein(:,5)==121,:) = [];
protein(protein(:,4) > 22,:) = [];

protein_water_select = [];
ias = 1;
for i = 1:size(protein_water_all,1)
    if min(dist(protein,protein_water_all(i,:)))< 2.9
        protein_water_select(ias,:) = protein_water_all(i,:);
        ias = ias + 1;
    end
end
        
protein = [protein
    protein_water_select];












