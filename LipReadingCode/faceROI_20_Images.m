clear all;
tic;

fid = fopen('fail.txt', 'wt');  % �إߤ@��.txt��(�˵L�k���L�B�ϰ쪺�ɮצW��)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};

for d = 1:6  % Ū�����O��Ƨ�
    
    cmddir = cmd{d};
    
    for user = 1 : 90   % Ū���ϥΪ̸�Ƨ�
        
        userdir = num2str(user);
        datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI'];  % Ū���W�٬�'1-90'����Ƨ�
%         mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI']);
        mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI_20_Images']);
%         mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
        datadir_20 = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\',];
        input_dir = dir(fullfile(datadir, '*.jpg'));  % �N���w��Ƨ��������w���ɦW���Ҧ��ɮ׸�ơA��Jinput_dir�}�C
        [x, y ] = size(input_dir);
        
        avarage = x/20;
        for n = 1 : 20
            [d,user,n]
            
            avarage_count = round(avarage * n);
            faceCatch20 = imread(fullfile(datadir, input_dir(avarage_count).name));
            faceCatch20 = imresize(faceCatch20,[96 96]);
            imwrite(faceCatch20,[datadir_20 'faceROI_20_Images\' input_dir(n).name]);
            total_image{(d-1)*90+user,n} = faceCatch20; 
        end
    end
end        