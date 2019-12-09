%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  data_layer_output(data_mat,lyr,output_dir,code_dir)

cd(output_dir);

%script1=[strcat('E:\OSDP\all\',list,'\',list,'_',num2str(asdf),'.txt')];
%save data.txt data_mat -ascii
eval(strcat('save data',num2str(lyr),'.txt data_mat -ascii'))
%save('data.txt', data_mat, '-ASCII');
cd(code_dir);


