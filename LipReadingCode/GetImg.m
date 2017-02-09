clc
clear

cmd = 'Walk';
inputdir = ['MatlabDataBase/six_cmd/videos/' cmd '/']; %影片路徑資料夾
outputd = ['MatlabDataBase/six_cmd/images/' cmd '/']; %圖片輸出資料夾
input_dir = dir(fullfile(inputdir, '*.mp4')); %讀取路淨資料夾下所有mp4檔的詳細資訊
[x, y] = size(input_dir); %x為檔案個數
for n = 1 : 1
n
    videonum = num2str(n); % 當前迴圈數轉字串
    outputdir = [outputd videonum '/'];
    mkdir(outputdir); % 新增一個輸出資料夾
    obj = VideoReader(fullfile(inputdir, input_dir(n).name)); %圖取當前影片詳細資訊
    objn = obj.NumberOfFrames; % 讀取當前影片Frame數
    objw = obj.Width; %影片寬
    objh = obj.Height; %影片高

    mov(1:objn) = struct('cdata',zeros(objw,objh,3,'uint8'),'colormap',[]); %建構cdata存720*960*3 Unit8圖片 ,colormap存空矩陣
    a = size(input_dir(n).name); %取得影片檔案名稱，例如 : Walk(1).mp4
    d = a(2)-4; %剪去檔名4最後個字(.mp4)
    f = input_dir(n).name(1:d); %存取剩下檔名(Walk(1))
    for k=1:objn
        mov(k).cdata = read(obj,k);
        output = mov(k).cdata;  %存取當前圖片
        tempStr = strcat(num2str(k),'.jpg'); % 存取檔名，例如:1.jpg

        if k < 10
            tempStr = [outputdir f '00' tempStr]; % 將檔名小於10的補00，例如 : 001.jpg
        else
            tempStr = [outputdir f  '0' tempStr]; %將檔名大於10的補0，例如 : 010.jpg
        end
        imwrite(output,tempStr); %將當前檔案寫入指定路徑
    end
end