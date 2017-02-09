clear all;
tic;
% img = imread('1/user11009.jpg');
% 建立膚色及唇色樣本模型
lipcolordir = 'LipcolorSample';
skincolordir = 'SkincolorSample';
[firstlipem, firstskinem] = lipgmm(lipcolordir, skincolordir);

fid = fopen('fail.txt', 'wt');  % 建立一個.txt檔(裝無法找到嘴唇區域的檔案名稱)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};
for d = 1:1  % 讀取數字資料夾
    cmddir = cmd{d};
for user = 1 : 1   % 讀取使用者資料夾
    userdir = num2str(user);
    datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\mouthROI\'];  % 讀取名稱為'1-90'的資料夾內的嘴巴區域影像
    input_dir = dir(fullfile(datadir, '*.jpg'));  % 將指定資料夾內的指定副檔名的所有檔案資料，丟入input_dir陣列
    [x, y ] = size(input_dir);
    for n = 1 : 1
        [user,n]

        lipROI = imread(fullfile(datadir, input_dir(n).name)); % 讀取指定資料夾內的指定檔案     
        
        lipROIsize = size(lipROI); % 抓取嘴巴區域的長寬
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
            
            testliparray = [];
            testskinarray = [];
            
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
            lipem = gmdistribution.fit(testliparray,2);
            
            testskinarray = double(testskinarray);
            skinem = gmdistribution.fit(testskinarray,2);
            
            sumlip = size(testliparray);
            sumlipcolorarr = [sumlipcolorarr; sumlip(1)];
            
            if i >=3
                if sumlipcolorarr(i) - sumlipcolorarr(i-1) < 0.78*(sumlipcolorarr(i-1) - sumlipcolorarr(i-2))
                    break;
                end
            end
        end
        
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
        %             imwrite(rboundinglip,[lipdir '\result\' inputlip_dir(n).name]);
        imwrite(rboundinglip,[datadir '\result\' input_dir(n).name]);
        % 轉成灰階影像
        graybb = rgb2gray(rboundinglip);
        imwrite(graybb,[datadir '\originalgray\' input_dir(n).name]);
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