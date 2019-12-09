%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  data_layer_xy_output(data_mat,lyr1,lyr2,output_dir,code_dir)

cd(output_dir);

%script1=[strcat('E:\OSDP\all\',list,'\',list,'_',num2str(asdf),'.txt')];
%save data.txt data_mat -ascii
eval(strcat('save data_',num2str(lyr1),'_',num2str(lyr2),'.txt data_mat -ascii'))
%save('data.txt', data_mat, '-ASCII');
cd(code_dir);


