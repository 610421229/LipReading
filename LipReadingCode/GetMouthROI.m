clear all;
tic;
% 建立膚色及唇色樣本模型
lipcolordir = 'LipcolorSample';
skincolordir = 'SkincolorSample';
[firstlipem, firstskinem] = lipgmm(lipcolordir, skincolordir);

fid = fopen('fail.txt', 'wt');  % 建立一個.txt檔(裝無法找到嘴唇區域的檔案名稱)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};
for d = 1:1  % 讀取指令資料夾
    cmddir = cmd{d};
for user = 1 : 90   % 讀取使用者資料夾
    userdir = num2str(user);
    datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\'];  % 讀取名稱為'1-90'的資料夾
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI']);
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','lipROI']);
    input_dir = dir(fullfile(datadir, '*.jpg'));  % 將指定資料夾內的指定副檔名的所有檔案資料，丟入input_dir陣列
    [x, y ] = size(input_dir);
    for n = 1 : x
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
            fprintf(fid,'%s\n',input_dir(n).name);
            continue; % 跳下一張圖片
        elseif fboxesx == 1 % 抓出臉部區域的點跟寬高
            left = fboxes(1);  % 左上角的 x 座標
            top = fboxes(2) - 10;   % 左上角的 y 座標 -10:怕切掉下嘴唇部分
            width = fboxes(3);  % 區域的寬
            height = fboxes(4) + 30;  % 區域的高 +30:怕切掉下嘴唇部分
        else % 找到不只一張臉時，取臉部面積最大的
            for i = fboxesx*2+1:fboxesx*3 % 從第一個的寬的欄位開始
                farea = fboxes(i)*fboxes(i + fboxesx); % 面積=寬*高
                if farea > maxarea % 目前臉部區域大於最大臉部區域，就取代並修改座標與最大面積值
                    maxarea = farea;
                    % 取臉部區域的座標點跟寬高
                    left = fboxes(i - 2*fboxesx);
                    top = fboxes(i - fboxesx) - 10;
                    width = fboxes(i);
                    height = fboxes(i + fboxesx) + 30;
                end
            end
        end
        faceROI = imcrop(img, [left top width height]);  % 將人臉區域剪下
        imwrite(faceROI,[datadir 'faceROI\' input_dir(n).name]);  % 將臉部區域存成檔案
        
        mouthDetector = vision.CascadeObjectDetector('Mouth');  % 嘴唇辨識器
        mboxes = step(mouthDetector, faceROI);  % 找出嘴唇
        lipROI = insertObjectAnnotation(faceROI, 'rectangle', mboxes, 'Mouth');
        [mboxesx, mboxesy] = size(mboxes);
        mouthleft = 0;
        if mboxesx ==0 % 找不到嘴巴就紀錄檔案名稱
            fprintf(fid,'%s\n',input_dir(n).name);
            continue; % 跳下一張圖
        end
        maxarea =2000;
        maxtop = 0; % 找最大的Y
        re = 0;   % 判斷是否抓取區域失敗的參數
        if mboxesx == 1
            mouthleft = mboxes(1);
            mouthtop = mboxes(2) -10;
            mouthwidth = mboxes(3);
            mouthheight = mboxes(4) +15;
            re = 1; % 代表抓取成功
        else
            for i = mboxesx+1:mboxesx*2  % 取Y欄位做為判斷標準
                larea = mboxes(i + mboxesx)*mboxes(i + mboxesx*2);  % 算嘴巴區域面積
                mouthright = mboxes(i - mboxesx) + mboxes(i + mboxesx);  % 找出最大的X
                %         if mboxes(i) > 2*height/3 && mboxes(i - mboxesx) < width/2 && mouthright > width/2 && larea > maxarea
                %     嘴唇區域判斷條件 : 最下面的點要>臉高的6.3/10，左右兩邊的點要在臉部中線的左右兩邊，可能的嘴唇區域中面積最大的
                if mboxes(i) > 6.3*height/10 && mboxes(i - mboxesx) < width/2 && mouthright > width/2 && larea > maxarea && mboxes(i) > maxtop
                    maxarea = larea;
                    maxtop = mboxes(i);
                    mouthleft = mboxes(i - mboxesx);
                    mouthtop = mboxes(i) - 10;
                    mouthwidth = mboxes(i + mboxesx);
                    mouthheight = mboxes(i + mboxesx*2) +15;
                    if mouthtop + mouthheight > height % 若是嘴唇的最低點超出臉部，則以臉部最低點為嘴唇最低點
                        mouthheight = height - mouthtop;
                    end
                    re = 1;
                end
            end
        end
        if re == 0
            fprintf(fid,'%s\n',input_dir(n).name);
            continue;
        end
        %  剪下最後的嘴巴區域並存檔
        mouthroi = [mouthleft mouthtop mouthwidth mouthheight];
        lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);   % 剪下最後的嘴巴區域
        imwrite(lipROI,[datadir 'mouthROI\' input_dir(n).name]);  % 將嘴巴區域存成檔案
        %figure,imshow(lipROI);
        lipROIsize = size(lipROI);
        modle = zeros(lipROIsize(1), lipROIsize(2), 3, 'uint8');
        % 將嘴唇影像轉成YCbCr色彩
        lipcolor = rgb2ycbcr(lipROI);
        slipcolor = size(lipcolor);
        % MxN的影像重新排列成M*Nx1的陣列
        rlipcolor = reshape(lipcolor, slipcolor(1)*slipcolor(2),3);
        rlipcolor = double(rlipcolor);
        % 計算影像中每個像素屬於膚色及唇色的機率
        flipe = pdf(firstlipem,rlipcolor);
        fskine = pdf(firstskinem,rlipcolor);
        
        % 利用YCbCr找出嘴唇
        for p=1:slipcolor(1)*slipcolor(2)
            if (fskine(p)>flipe(p))
                rlipcolor(p,:) = [0.0627*255 0.5020*255 0.5020*255];
            end
        end
        rlipcolor = reshape(rlipcolor, slipcolor(1), slipcolor(2), 3);
        rlipcolor = uint8(rlipcolor);
        
        I = ycbcr2rgb(rlipcolor);  % 轉回RGB影像
        
        %figure, imshow(I);
        
        g = rgb2gray(I);  % RGB轉灰階圖
        bw = g > 0;  % 將灰階圖二值化
        
        % 二值圖做侵蝕
        se = strel('disk',4);
        erobw = imerode(bw,se);
        
        bcc = erobw;
        cc = bwconncomp(bcc);
        numpixels = cellfun(@length,cc.PixelIdxList);  % 找出影像中的物件(白色區域)的結果
        [numsort, idxsort] = sort(numpixels);
        [biggest, bidx] = max(numpixels);
        second = 0; sid = 0;
        len = length(numsort);
        %figure, imshow(bcc);
        if len >= 2
            second = numsort(len-1);
            sid = idxsort(len-1);
        end
        %             保留最大塊的白色區域
        for i = 1 : max(idxsort)
            if i == bidx
                bcc(cc.PixelIdxList{i}) = 1;
                %     elseif i == sid && second > 300
                %         bcc(cc.PixelIdxList{i}) = 1;
            else
                bcc(cc.PixelIdxList{i}) = 0;
            end
        end
        %figure, imshow(bcc);
        doublelip = double(lipROI);
        modle(:,:,1) = bcc .* doublelip(:,:,1);
        modle(:,:,2) = bcc .* doublelip(:,:,2);
        modle(:,:,3) = bcc .* doublelip(:,:,3);
        
        % 利用侵蝕後的結果找出新的唇色模型與非唇色模型建立
        liparray = [];
        skinarray = [];
        
        for i = 1:lipROIsize(1)    %mouthheight+1
            for j = 1:lipROIsize(2)   %mouthwidth+1
                if modle(i,j) ~= 0
                    liparray = [liparray;doublelip(i,j,1), doublelip(i,j,2), doublelip(i,j,3)];
                else
                    skinarray = [skinarray;doublelip(i,j,1), doublelip(i,j,2), doublelip(i,j,3)];
                end
            end
        end
        % [asx, asy] = size(liparray);
        % if asx == 0 || asx/(lipROIsize(1)*lipROIsize(2)) < 0.01
        %     return;
        % end
        liparray = double(liparray);
        lipem = gmdistribution.fit(liparray,1);
        
        skinarray = double(skinarray);
        skinem = gmdistribution.fit(skinarray,1);
        
        sumlipcolorarr = [];
        % 迴圈色唇色偵測
        for i = 1 : 10
            testlip = lipROI;
            stestlip = size(testlip);
            rtestlip = reshape(testlip, stestlip(1)*stestlip(2),3);
            rtestlip = double(rtestlip);
            
            % pdf機率模型 - 唇色
            lipe = pdf(lipem,rtestlip);
            % pdf機率模型 - 膚色
            skine = pdf(skinem,rtestlip);
            % 門檻值theta
            theta = [5e-6];
            
            % 唇色判定
            for p=1:stestlip(1)*stestlip(2)
                if (skine(p)>lipe(p))
                    rtestlip(p,:) = [0 0 0];
                end
            end
            
            result  = reshape(rtestlip, stestlip(1), stestlip(2), 3);
            result = uint8(result);
            
            testliparray = []; % 新的唇色樣本陣列
            testskinarray = []; % 新的膚色樣本陣列
            
            for j = 1:stestlip(1)
                for k = 1:stestlip(2)
                    if result(j,k) ~= 0
                        testliparray = [testliparray;doublelip(j,k,1), doublelip(j,k,2), doublelip(j,k,3)];
                    else
                        testskinarray = [testskinarray;doublelip(j,k,1), doublelip(j,k,2), doublelip(j,k,3)];
                    end
                end
            end
            
            testliparray = double(testliparray);
            lipem = gmdistribution.fit(testliparray,2); % 新的唇色模型
            
            testskinarray = double(testskinarray);
            skinem = gmdistribution.fit(testskinarray,2); % 新的膚色模型
            
            sumlip = size(testliparray);
            sumlipcolorarr = [sumlipcolorarr; sumlip(1)];
            
            if i >=3
                if sumlipcolorarr(i) - sumlipcolorarr(i-1) < 0.78*(sumlipcolorarr(i-1) - sumlipcolorarr(i-2))
                    break;
                end
            end
        end
        
        %figure, imshow(result);
        
        % 迴圈式結果侵蝕
        g = rgb2gray(result);  % RGB轉灰階圖
        bw = g > 0;  % 將灰階圖二值化
        
        % 二值圖做closing
        se = strel('disk',1);
        erobw = imclose(bw,se);
        
        bcc = erobw;
        cc = bwconncomp(bcc);
        numpixels = cellfun(@length,cc.PixelIdxList);
        [numsort, idxsort] = sort(numpixels);
        [biggest, bidx] = max(numpixels);
        second = 0; sid = 0;
        len = length(numsort);
        
        if len >= 2
            second = numsort(len-1);
            sid = idxsort(len-1);
        end
        
        for i = 1 : max(idxsort)
            if i == bidx
                bcc(cc.PixelIdxList{i}) = 1;
                %     elseif i == sid && second > 300
                %         bcc(cc.PixelIdxList{i}) = 1;
            else
                bcc(cc.PixelIdxList{i}) = 0;
            end
        end
        
        %figure, imshow(bcc);
        
        % 找出BoundingBox的區域
        topwp = 0; t = 0;
        underwp = 0; u = 0;
        leftwp = 0; l = 0;
        rightwp = 0; r = 0;
        % 找出最上面的點
        for i = 1 : lipROIsize(1)   %mouthheight+1
            for j = 1 : lipROIsize(2)   %mouthwidth+1
                if bcc(i,j) ~=0 && t == 0
                    topwp = i;
                    topx = j;
                    t = 1;
                elseif t == 1
                    break;
                end
            end
        end
        % 找出最下面的點
        for i = lipROIsize(1) : -1 : 1   %mouthheight+1 : -1 : 1
            for j = lipROIsize(2) : -1 : 1   %mouthwidth+1 : -1 : 1
                if bcc(i,j) ~=0 && u == 0
                    underwp = i;
                    u = 1;
                elseif u == 1
                    break;
                end
            end
        end
        % 找出最左邊的點
        for j = 1 : lipROIsize(2)   %mouthwidth+1
            for i = 1 : lipROIsize(1)   %mouthheight+1
                if bcc(i,j) ~=0 && l == 0
                    leftwp = j;
                    l = 1;
                elseif l == 1
                    break;
                end
            end
        end
        % 找出最右邊的點
        for j = lipROIsize(2) : -1 : 1   %mouthwidth+1 : -1 : 1
            for i = lipROIsize(1) : -1 : 1   %mouthheight+1 : -1 : 1
                if bcc(i,j) ~=0 && r == 0
                    rightwp = j;
                    r = 1;
                elseif r == 1
                    break;
                end
            end
        end
        
        wpwidth = rightwp - leftwp;
        wpheight = underwp - topwp;
        
        boundinglip = imcrop(lipROI, [leftwp topwp wpwidth wpheight]);
        
        %figure, imshow(boundinglip);
        
        boundinglipgray = rgb2gray(boundinglip);
        [bbsize1, bbsize2] = size(bcc);
        % bbsize1 = bbsize(1);
        % bbsize2 = bbsize(2);
        edgeofROI = bwtraceboundary(bcc,[topwp, topx],'W',8);
        % cannyedgeofROI = edge(bcc, 'canny');
        bwtsize = size(edgeofROI);
        bwtrace = zeros(bbsize1, bbsize2, 1);
        for i = 1 : bwtsize(1)
            bwtrace(edgeofROI(i),edgeofROI(i+bwtsize(1))) = 1;
        end
        
        %figure, imshow(bwtrace);
        
        edgepointsxy = [];
        upskedgey = [];
        upskedgex = [];
        lowskedgey = [];
        lowskedgex = [];
        % 對boundingbox內的唇做二次曲線擬合，y = ax^2+bx+c，先判斷像素是否為唇色，若屬於唇色則留下
        for i = 1 : bbsize1
            for j = 1 : bbsize2
                points = [lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
                points = double(points);
                % 計算BoundingBox內的像素屬於唇色的機率
                flipe = pdf(lipem, points);
                fskine = pdf(skinem, points);
                if bwtrace(i,j) == 1 && i <= bbsize1/2 % 上唇
                    if fskine < flipe
                        %                 upskedge = [upskedge; lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
                        upskedgey = [upskedgey; i];                           % y = [y1;y2;...;yn]
                        upskedgex = [upskedgex; j^2, j, 1];              % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
                    else
                        bwtrace(i,j) = 0;
                    end
                elseif bwtrace(i,j) == 1 && i > bbsize1/2 % 下唇
                    if fskine < flipe
                        %                 lowskedge = [lowskedge; lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
                        lowskedgey = [lowskedgey; i];                           % y = [y1;y2;...;yn]
                        lowskedgex = [lowskedgex; j^2, j, 1];              % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
                    else
                        bwtrace(i,j) = 0;
                    end
                end
            end
        end
        
        % 判斷屬於唇色的點是否為曲線上的點-上唇
        f = 0;
        upr = 0;
        c = 3;
        updelta = 3;
        upxsize = size(upskedgex);    % Eskin的數量
        uprx = upskedgex;             % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
        upry = upskedgey;             % y = [y1;y2;...;yn]
        for t = 1 : 10
            uprxsize = size(uprx);
            upx = pinv(uprx);
            upu = upx*upry;
            if uprxsize <= upr
                f = f + 1;
                if f > c
                    break;
                end
            else
                f = 0;
            end
            Rx = [];
            Ry = [];
            for i = 1 : upxsize(1)
                p = abs(upu(1)*upskedgex(i) + upu(2)*upskedgex(i+upxsize(1)) + upu(3) - upskedgey(i));
                if p < updelta
                    Rx = [Rx; upskedgex(i), upskedgex(i+upxsize(1)), upskedgex(i+upxsize(1)*2)];
                    Ry = [Ry; upskedgey(i)];
                end
            end
            uprx = Rx;
            upry = Ry;
            upr = uprxsize(1);
        end
        
        upr = size(uprx);
        uplip = zeros(bbsize1, bbsize2, 1);
        for i = 1:upr(1)
            uplip(upry(i), uprx(i+upr(1))) = 1;
        end
        uplip = double(uplip);
        
        %figure, imshow(uplip);
        
        % 判斷屬於唇色的點是否為曲線上的點-下唇
        f = 0;
        lowdelta = 1;
        lowr = 0;
        lowxsize = size(lowskedgex);
        lowrx = lowskedgex;
        lowry = lowskedgey;
        for t = 1 : 10
            lowrxsize = size(lowrx);
            lowx = pinv(lowrx);
            lowu = lowx*lowry;
            if lowrxsize <= lowr
                f = f + 1;
                if f > c
                    break;
                end
            else
                f = 0;
            end
            Rx = [];
            Ry = [];
            for i = 1 : lowxsize(1)
                p = abs(lowu(1)*lowskedgex(i) + lowu(2)*lowskedgex(i+lowxsize(1)) + lowu(3) - lowskedgey(i));
                if p < lowdelta
                    Rx = [Rx; lowskedgex(i), lowskedgex(i+lowxsize(1)), lowskedgex(i+lowxsize(1)*2)];
                    Ry = [Ry; lowskedgey(i)];
                end
            end
            lowrx = Rx;
            lowry = Ry;
            lowr = lowrxsize(1);
        end
        lowr = size(lowrx);
        lowlip = zeros(bbsize1, bbsize2, 1);
        for i = 1:lowr(1)
            lowlip(lowry(i), lowrx(i+lowr(1))) = 1;
        end
        lowlip = double(lowlip);
        
        %figure, imshow(lowlip);
        
        edgepointsxy = [];
        for i = 1 : upr(1) + lowr(1)
            if i <= upr(1)
                edgepointsxy = [edgepointsxy; uprx(i+upr(1)), upry(i)];
            elseif i >upr(1)
                edgepointsxy = [edgepointsxy; lowrx(i-upr(1)+lowr(1)), lowry(i-upr(1))];
            end
        end
        sizeofp = size(edgepointsxy);
        tr = zeros(bbsize1, bbsize2, 1);
        for i = 1:sizeofp(1)
            tr(edgepointsxy(i+sizeofp(1)),edgepointsxy(i)) = 1;
        end
        tr = logical(tr);
        
        % 找出BoundingBox的區域
        topwp = 0; t = 0;
        underwp = 0; u = 0;
        leftwp = 0; l = 0;
        rightwp = 0; r = 0;
        % find the top point
        for i = 1 : bbsize1   %mouthheight+1
            for j = 1 : bbsize2   %mouthwidth+1
                if tr(i,j) ~=0 && t == 0
                    topwp = i;
                    t = 1;
                elseif t == 1
                    break;
                end
            end
        end
        % find the under point
        for i = bbsize1 : -1 : 1   %mouthheight+1 : -1 : 1
            for j = bbsize2 : -1 : 1   %mouthwidth+1 : -1 : 1
                if tr(i,j) ~=0 && u == 0
                    underwp = i;
                    u = 1;
                elseif u == 1
                    break;
                end
            end
        end
        % find the left point
        for j = 1 : bbsize2   %mouthwidth+1
            for i = 1 : bbsize1   %mouthheight+1
                if tr(i,j) ~=0 && l == 0
                    leftwp = j;
                    l = 1;
                elseif l == 1
                    break;
                end
            end
        end
        % find the right point
        for j =bbsize2 : -1 : 1   %mouthwidth+1 : -1 : 1
            for i = bbsize1 : -1 : 1   %mouthheight+1 : -1 : 1
                if tr(i,j) ~=0 && r == 0
                    rightwp = j;
                    r = 1;
                elseif r == 1
                    break;
                end
            end
        end
        
        wpwidth = rightwp - leftwp;
        wpheight = underwp - topwp;
        % 剪出最後的嘴唇區域
        rboundinglip = imcrop(lipROI, [leftwp topwp wpwidth wpheight]);
        imwrite(rboundinglip,[datadir 'lipROI\' input_dir(n).name]);
        %figure, imshow(rboundinglip);
%         %             imwrite(rboundinglip,[lipdir '\result\' inputlip_dir(n).name]);
%         imwrite(rboundinglip,[datadir '\result\' input_dir(n).name]);
        % 轉成灰階影像
        graybb = rgb2gray(rboundinglip);
%         imwrite(graybb,[datadir '\originalgray\' input_dir(n).name]);
        %
        %             grayimg = 'newdatabase/data/1/ok/originalgray/user261001.jpg';
        %             grayimg = imread(grayimg);
        %
        %             xgray = double(grayimg);
        %             mx = mean(xgray(:));
        %             sx = std(xgray(:));
        %             ygray = double(graybb);
        %             my = mean(ygray(:));    % 算平均值
        %             sy = std(ygray(:));     % 算標準差
        %
        %             normaly = (ygray-my)/sy*sx+mx;
        %             normaly = uint8(normaly);
        %             imwrite(normaly,[lipdir '\luminanceNor\' inputlip_dir(n).name]);
        %             imwrite(normaly,['newdatabase\testdata\luminanceNor\' userdir '\' facedir '\' inputlip_dir(n).name]);
    end
    
    % imwrite(boundinglip, '1/user11008.jpg');
end
end
fclose(fid);
toc;