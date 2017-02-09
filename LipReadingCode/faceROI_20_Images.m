clear all;
tic;

fid = fopen('fail.txt', 'wt');  % 建立一個.txt檔(裝無法找到嘴唇區域的檔案名稱)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};

for d = 1:6  % 讀取指令資料夾
    
    cmddir = cmd{d};
    
    for user = 1 : 90   % 讀取使用者資料夾
        
        userdir = num2str(user);
        datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI'];  % 讀取名稱為'1-90'的資料夾
%         mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI']);
        mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI_20_Images']);
%         mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
        datadir_20 = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\',];
        input_dir = dir(fullfile(datadir, '*.jpg'));  % 將指定資料夾內的指定副檔名的所有檔案資料，丟入input_dir陣列
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