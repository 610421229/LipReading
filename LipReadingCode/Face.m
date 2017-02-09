clear all;
tic;

fid = fopen('fail.txt', 'wt');  % 建立一個.txt檔(裝無法找到嘴唇區域的檔案名稱)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};

for d = 2:6  % 讀取指令資料夾
    
    cmddir = cmd{d};
    
    for user = 1 : 90   % 讀取使用者資料夾
        
        userdir = num2str(user);
        datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\'];  % 讀取名稱為'1-90'的資料夾
        mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI']);
        mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
        input_dir = dir(fullfile(datadir, '*.jpg'));  % 將指定資料夾內的指定副檔名的所有檔案資料，丟入input_dir陣列
        [x, y ] = size(input_dir);
        
        for n = 1 : x
                
            [user,n]
            img = imread(fullfile(datadir, input_dir(n).name)); % 讀取指定資料夾內的指定檔案
            img = imresize(img,[360 480]);
            faceDetector =vision.CascadeObjectDetector;  % 人臉辨識器
            fboxes = step(faceDetector, img);  % 找出人臉
            %  fboxes[a,b,c,d]: a->左上角的x座標 b->左上角的y座標 c->臉部區域的寬 d->臉部區域的高
            faceROI = insertObjectAnnotation(img, 'rectangle', fboxes, 'Face');
            [fboxesx, fboxesy] = size(fboxes);
            farea = 0; % 臉部面積參數
            maxarea = 0; % 最大臉部面積參數
            
                if fboxesx == 0 % 找不到人臉就紀錄檔案名稱
                    %fprintf(fid,'%s\n',input_dir(n).name);
                    %continue; % 跳下一張圖片
                
                elseif fboxesx == 1 % 抓出臉部區域的點跟寬高
                    
                    left = fboxes(1);  % 左上角的 x 座標
                    top = fboxes(2) - 10;   % 左上角的 y 座標 -10:怕切掉下嘴唇部分
                    width = fboxes(3);  % 區域的寬
                    height = fboxes(4) + 30;  % 區域的高 +30:怕切掉下嘴唇部分
                    
                else % 找到不只一張臉時，取臉部面積最大的
                    
                    for i = fboxesx*2+1 : fboxesx*3 % 從第一個的寬的欄位開始
                        
                            farea = fboxes(i)*fboxes(i + fboxesx); % 面積=寬*高
                        
                            if farea > maxarea % 目前臉部區域大於最大臉部區域，就取代並修改座標與最大面積值
                            
                            maxarea = farea;
                            % 取臉部區域的座標點跟寬高
                            left = fboxes(i - 2*fboxesx);
                            top = fboxes(i - fboxesx) - 10;
                            width = fboxes(i);
                            height = fboxes(i + fboxesx) ;
                            
                            end
                    end
                end
        faceROI = imcrop(img, [left top+height*2/3 width height/3]);% 將人臉區域剪下
        %faceROI = imcrop(faceROI, [left top+height*2/3 width height]);
        imwrite(faceROI,[datadir 'faceROI\' input_dir(n).name]);  % 將臉部區域存成檔案
        
        mouthDetector = vision.CascadeObjectDetector('Mouth');  % 嘴唇辨識器
        mboxes = step(mouthDetector, faceROI);  % 找出嘴唇
        lipROI = insertObjectAnnotation(faceROI, 'rectangle', mboxes, 'Mouth');
        [mboxesx, mboxesy] = size(mboxes);
        mouthleft = 0;
        maxarea =2000;
        maxtop = 0; % 找最大的Y
        re = 0;   % 判斷是否抓取區域失敗的參數
        
        if mboxesx ==0 % 找不到嘴巴就紀錄檔案名稱
            lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);
            %fprintf(fid,'%s\n',input_dir(n).name);
            %continue; % 跳下一張圖
            
        elseif mboxesx == 1
            mouthleft = mboxes(1);
            mouthtop = mboxes(2) -10;
            mouthwidth = mboxes(3);
            mouthheight = mboxes(4) +15;
            re = 1; % 代表抓取成功
        else

            for i = mboxesx+1:mboxesx*2  % 取Y欄位做為判斷標準
                larea = mboxes(i + mboxesx)*mboxes(i + mboxesx*2);  % 算嘴巴區域面積
            
                if larea>maxarea;
                    maxarea = larea;
                    mouthleft = mboxes(i - mboxesx);
                    mouthtop = mboxes(i) - 10;
                    mouthwidth = mboxes(i + mboxesx);
                    mouthheight = mboxes(i + mboxesx*2) ;
                end
            end
        end

        mouthroi = [mouthleft mouthtop mouthwidth mouthheight];
        lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);   % 剪下最後的嘴巴區域
        imwrite(lipROI,[datadir 'mouthROI\' input_dir(n).name]);  % 將嘴巴區域存成檔案

        end
    end
end
fclose(fid);
toc;