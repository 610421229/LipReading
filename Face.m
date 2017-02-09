clear all;
tic;
%% 建立膚色及唇色樣本模型
% lipcolordir = 'LipcolorSample';
% skincolordir = 'SkincolorSample';
% [firstlipem, firstskinem] = lipgmm(lipcolordir, skincolordir);
%%
fid = fopen('fail.txt', 'wt');  % 建立一個.txt檔(裝無法找到嘴唇區域的檔案名稱)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};
for d = 1:1  % 讀取指令資料夾
    cmddir = cmd{d};
for user = 1 : 1   % 讀取使用者資料夾
    userdir = num2str(user);
    datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\'];  % 讀取名稱為'1-90'的資料夾
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','lipROI']);
    input_dir = dir(fullfile(datadir, '*.jpg'));  % 將指定資料夾內的指定副檔名的所有檔案資料，丟入input_dir陣列
    [x, y ] = size(input_dir);
     for n = 28 : 29
        [user,n]
        img = imread(fullfile(datadir, input_dir(n).name)); % 讀取指定資料夾內的指定檔案
        faceDetector =vision.CascadeObjectDetector;  % 人臉辨識器
        fboxes = step(faceDetector, img);  % 找出人臉
        %  fboxes[a,b,c,d]: a->左上角的x座標 b->左上角的y座標 c->臉部區域的寬 d->臉部區域的高
        faceROI = insertObjectAnnotation(img, 'rectangle', fboxes, 'Face');
        [fboxesx, fboxesy] = size(fboxes);
        farea = 0; % 臉部面積參數
        maxarea = 0; % 最大臉部面積參數
        if fboxesx == 0 % 找不到人臉就紀錄檔案名稱