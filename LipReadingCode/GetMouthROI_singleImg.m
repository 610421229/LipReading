close all
tic;
filename = 'Drink (01)001.jpg';
dir = ['MatlabDataBase\six_cmd\images\Drink\1\' filename];
roidir = ['MatlabDataBase\six_cmd\images\Drink\1\mouthROI\' filename];
% resultdir = ['newdatabase\testdata\user135\4\result\' filename];
% graydir = ['newdatabase\testdata\user135\4\originalgray\' filename];
% nordir = ['newdatabase\testdata\user135\4\luminanceNor\' filename];
% lumdir = ['newdatabase\testdata\luminanceNor\user135\4\' filename];
img = imread(dir);
% path = '3/used/user23006.jpg';
% rpath = '3/result/user230062.jpg';
% img = imread(path);
faceDetector = vision.CascadeObjectDetector;  % 人臉辨識器
fboxes = step(faceDetector, img);  % 找出人臉
faceROI = insertObjectAnnotation(img, 'rectangle', fboxes, 'Face');
figure, imshow(faceROI);
[fboxesx, fboxesy] = size(fboxes);
farea = 0;
maxarea = 0;
if fboxesx == 1
    left = fboxes(1);  % 左上角的 x 座標
    top = fboxes(2) - 10;   % 左上角的 y 座標
    width = fboxes(3);  % 區域的寬
    height = fboxes(4) + 30;  % 區域的高 +30:怕切掉下嘴唇部分
else
    for i = fboxesx*2+1:fboxesx*3
        farea = fboxes(i)*fboxes(i + fboxesx);
        if farea > maxarea
            maxarea = farea;
            left = fboxes(i - 2*fboxesx);
            top = fboxes(i - fboxesx) - 10;
            width = fboxes(i);
            height = fboxes(i + fboxesx) + 30;
        end
    end
end
faceROI = imcrop(img, [left top width height]);  % 將人臉區域剪下
% figure, imshow(faceROI);

mouthDetector = vision.CascadeObjectDetector('Mouth');  % 嘴唇辨識器
mboxes = step(mouthDetector, faceROI);  % 找出嘴唇
lipROI = insertObjectAnnotation(faceROI, 'rectangle', mboxes, 'Mouth');
figure, imshow(lipROI);
[mboxesx, mboxesy] = size(mboxes);
mouthleft = 0;
if mboxesx ==0
    return;
end
maxarea = 5000;
maxtop = 0;
if mboxesx == 1
    mouthleft = mboxes(1) - 20;
    mouthtop = mboxes(2) - 20;
    mouthwidth = mboxes(3) + 40;
    mouthheight = mboxes(4) + 45;
else
    for i = mboxesx+1:mboxesx*2
        larea = mboxes(i + mboxesx)*mboxes(i + mboxesx*2);
        mouthright = mboxes(i - mboxesx) + mboxes(i + mboxesx);
        %         if mboxes(i) > 2*height/3 && mboxes(i - mboxesx) < width/2 && mouthright > width/2 && larea > maxarea
        if mboxes(i) > 6.3*height/10 && mboxes(i - mboxesx) < width/2 && mouthright > width/2 && larea > maxarea && mboxes(i) > maxtop
            maxarea = larea;
            maxtop = mboxes(i);
            mouthleft = mboxes(i - mboxesx) - 20;
            mouthtop = mboxes(i) - 20;
            mouthwidth = mboxes(i + mboxesx) + 40;
            mouthheight = mboxes(i + mboxesx*2) + 45;
            if mouthtop + mouthheight > height
                mouthheight = height - mouthtop;
            end
        end
    end
end
% if re == 0
%     return;
% end
mouthroi = [mouthleft mouthtop mouthwidth mouthheight];
lipROIcheck = insertObjectAnnotation(faceROI, 'rectangle', mouthroi, 'Mouth');
figure, imshow(lipROIcheck);
lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);
imwrite(lipROI,roidir);
% % figure, imshow(lipROI);
% % figure, subplot(2, 3, 1), imshow(lipROI);
% figure, subplot(1, 2, 1), imshow(lipROI);
% 
% lipROIsize = size(lipROI);
% modle = zeros(lipROIsize(1), lipROIsize(2), 3, 'uint8');
% 
% lipcolor = rgb2ycbcr(lipROI);
% slipcolor = size(lipcolor);
% rlipcolor = reshape(lipcolor, slipcolor(1)*slipcolor(2),3);
% rlipcolor = double(rlipcolor);
% 
% lipcolordir = 'lipcolor';
% skincolordir = 'skincolor';
% [firstlipem, firstskinem] = lipgmm(lipcolordir, skincolordir);
% flipe = pdf(firstlipem,rlipcolor);
% fskine = pdf(firstskinem,rlipcolor);
% 
% % 利用YCbCr找出嘴唇
% for p=1:slipcolor(1)*slipcolor(2)
%     if (fskine(p)>flipe(p))
%         rlipcolor(p,:) = [0.0627*255 0.5020*255 0.5020*255];
%     end
% end
% rlipcolor = reshape(rlipcolor, slipcolor(1), slipcolor(2), 3);
% rlipcolor = uint8(rlipcolor);
% 
% I = ycbcr2rgb(rlipcolor);
% g = rgb2gray(I);  % RGB轉灰階圖
% bw = g > 0;  % 將灰階圖二值化
% 
% % 二值圖做侵蝕
% se = strel('disk',4);
% erobw = imerode(bw,se);
% 
% bcc = erobw;
% cc = bwconncomp(bcc);
% numpixels = cellfun(@length,cc.PixelIdxList);
% [numsort, idxsort] = sort(numpixels);
% [biggest, bidx] = max(numpixels);
% second = 0; sid = 0;
% len = length(numsort);
% 
% if len >= 2
%     second = numsort(len-1);
%     sid = idxsort(len-1);
% end
% 
% for i = 1 : max(idxsort)
%     if i == bidx
%         bcc(cc.PixelIdxList{i}) = 1;
% %     elseif i == sid && second > 300
% %         bcc(cc.PixelIdxList{i}) = 1;
%     else
%         bcc(cc.PixelIdxList{i}) = 0;
%     end
% end
% 
% figure, imshow(bcc);
% 
% doublelip = double(lipROI);
% modle(:,:,1) = bcc .* doublelip(:,:,1);
% modle(:,:,2) = bcc .* doublelip(:,:,2);
% modle(:,:,3) = bcc .* doublelip(:,:,3);
% % figure, imshow(modle);
% 
% % 唇色模型與非唇色模型建立
% liparray = [];
% skinarray = [];
% 
% for i = 1:lipROIsize(1)    %mouthheight+1
%     for j = 1:lipROIsize(2)   %mouthwidth+1
%         if modle(i,j) ~= 0
%             liparray = [liparray;doublelip(i,j,1), doublelip(i,j,2), doublelip(i,j,3)];
%         else
%             skinarray = [skinarray;doublelip(i,j,1), doublelip(i,j,2), doublelip(i,j,3)];
%         end
%     end
% end
% % [asx, asy] = size(liparray);
% % if asx == 0 || asx/(lipROIsize(1)*lipROIsize(2)) < 0.01
% %     return;
% % end
% liparray = double(liparray);
% lipem = gmdistribution.fit(liparray,1);
% 
% skinarray = double(skinarray);
% skinem = gmdistribution.fit(skinarray,1);
% 
% sumlipcolorarr = [];
% 
% for i = 1 : 10
%     testlip = lipROI;
%     stestlip = size(testlip);
%     rtestlip = reshape(testlip, stestlip(1)*stestlip(2),3);
%     rtestlip = double(rtestlip);
% 
%     % pdf機率模型 - 唇色
%     lipe = pdf(lipem,rtestlip);
%     % pdf機率模型 - 膚色
%     skine = pdf(skinem,rtestlip);
%     % 門檻值theta
%     theta = [5e-6];
% 
%     % 唇色判定
%     for p=1:stestlip(1)*stestlip(2)
%         if (skine(p)>lipe(p))
%             rtestlip(p,:) = [0 0 0];
%         end
%     end
% 
%     result  = reshape(rtestlip, stestlip(1), stestlip(2), 3);
%     result = uint8(result);
%     
%     testliparray = [];
%     testskinarray = [];
%     
%     for j = 1:stestlip(1)
%         for k = 1:stestlip(2)
%             if result(j,k) ~= 0
%                 testliparray = [testliparray;doublelip(j,k,1), doublelip(j,k,2), doublelip(j,k,3)];
%             else
%                 testskinarray = [testskinarray;doublelip(j,k,1), doublelip(j,k,2), doublelip(j,k,3)];
%             end
%         end
%     end
%     
%     testliparray = double(testliparray);
%     lipem = gmdistribution.fit(testliparray,1);
% 
%     testskinarray = double(testskinarray);
%     skinem = gmdistribution.fit(testskinarray,1);
%     
%     sumlip = size(testliparray);
%     sumlipcolorarr = [sumlipcolorarr; sumlip(1)];
%     
%     if i >=3
%         if sumlipcolorarr(i) - sumlipcolorarr(i-1) < 0.83*(sumlipcolorarr(i-1) - sumlipcolorarr(i-2))
%             break;
%         end
%     end
% end
% figure, imshow(result);
% 
% % 迴圈式結果侵蝕
% g = rgb2gray(result);  % RGB轉灰階圖
% bw = g > 0;  % 將灰階圖二值化
% 
% % 二值圖做closing
% se = strel('disk',1);
% erobw = imclose(bw,se);
% 
% bcc = erobw;
% cc = bwconncomp(bcc);
% numpixels = cellfun(@length,cc.PixelIdxList);
% [numsort, idxsort] = sort(numpixels);
% [biggest, bidx] = max(numpixels);
% second = 0; sid = 0;
% len = length(numsort);
% 
% if len >= 2
%     second = numsort(len-1);
%     sid = idxsort(len-1);
% end
% 
% for i = 1 : max(idxsort)
%     if i == bidx
%         bcc(cc.PixelIdxList{i}) = 1;
% %     elseif i == sid && second > 300
% %         bcc(cc.PixelIdxList{i}) = 1;
%     else
%         bcc(cc.PixelIdxList{i}) = 0;
%     end
% end
% % % figure, imshow(result);
% % subplot(1, 2, 2), imshow(bcc);
% 
% % 找出BoundingBox的區域
% topwp = 0; t = 0;
% underwp = 0; u = 0;
% leftwp = 0; l = 0;
% rightwp = 0; r = 0;
% 
% for i = 1 : lipROIsize(1)   %mouthheight+1
%     for j = 1 : lipROIsize(2)   %mouthwidth+1
%         if bcc(i,j) ~=0 && t == 0
%             topwp = i;
%             topx = j;
%             t = 1;
%         elseif t == 1
%             break;
%         end
%     end
% end
% 
% for i = lipROIsize(1) : -1 : 1   %mouthheight+1 : -1 : 1
%     for j = lipROIsize(2) : -1 : 1   %mouthwidth+1 : -1 : 1
%         if bcc(i,j) ~=0 && u == 0
%             underwp = i;
%             u = 1;
%         elseif u == 1
%             break;
%         end
%     end
% end
% 
% for j = 1 : lipROIsize(2)   %mouthwidth+1
%     for i = 1 : lipROIsize(1)   %mouthheight+1
%         if bcc(i,j) ~=0 && l == 0
%             leftwp = j;
%             l = 1;
%         elseif l == 1
%             break;
%         end
%     end
% end
% 
% for j = lipROIsize(2) : -1 : 1   %mouthwidth+1 : -1 : 1
%     for i = lipROIsize(1) : -1 : 1   %mouthheight+1 : -1 : 1
%         if bcc(i,j) ~=0 && r == 0
%             rightwp = j;
%             r = 1;
%         elseif r == 1
%             break;
%         end
%     end
% end
% 
% wpwidth = rightwp - leftwp;
% wpheight = underwp - topwp;
% 
% boundinglip = imcrop(lipROI, [leftwp topwp wpwidth wpheight]);
% boundinglipgray = rgb2gray(boundinglip);
% figure, imshow(boundinglip);
% [bbsize1, bbsize2] = size(bcc);
% % bbsize1 = bbsize(1);
% % bbsize2 = bbsize(2);
% edgeofROI = bwtraceboundary(bcc,[topwp, topx],'W',8);
% % cannyedgeofROI = edge(bcc, 'canny');
% bwtsize = size(edgeofROI);
% bwtrace = zeros(bbsize1, bbsize2, 1);
% for i = 1 : bwtsize(1)
%     bwtrace(edgeofROI(i),edgeofROI(i+bwtsize(1))) = 1;
% end
% % subplot(2, 3, 3), imshow(bwtrace);
% % figure, subplot(1, 3, 1), imshow(bwtrace);
% figure, imshow(bwtrace);
% edgepointsxy = [];
% % upskedge = [];
% upskedgey = [];
% upskedgex = [];
% % lowskedge = [];
% lowskedgey = [];
% lowskedgex = [];
% % y = ax^2+bx+c
% for i = 1 : bbsize1
%     for j = 1 : bbsize2
%         points = [lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
%         points = double(points);
%         flipe = pdf(lipem, points);
%         fskine = pdf(skinem, points);
%         if bwtrace(i,j) == 1 && i <= bbsize1/2
%             if fskine < flipe
% %                 upskedge = [upskedge; lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
%                 upskedgey = [upskedgey; i];                           % y = [y1;y2;...;yn]
%                 upskedgex = [upskedgex; j^2, j, 1];              % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
%             else
%                 bwtrace(i,j) = 0;
%             end
%         elseif bwtrace(i,j) == 1 && i > bbsize1/2
%             if fskine < flipe
% %                 lowskedge = [lowskedge; lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
%                 lowskedgey = [lowskedgey; i];                           % y = [y1;y2;...;yn]
%                 lowskedgex = [lowskedgex; j^2, j, 1];              % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
%             else
%                 bwtrace(i,j) = 0;
%             end
%         end
%     end
% end
% 
% f = 0;
% upr = 0;
% c = 3;
% updelta = 3;
% upxsize = size(upskedgex);    % Eskin的數量
% uprx = upskedgex;             % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
% upry = upskedgey;             % y = [y1;y2;...;yn]
% for t = 1 : 10
%     uprxsize = size(uprx);
%     upx = pinv(uprx);
%     upu = upx*upry;
%     if uprxsize <= upr
%         f = f + 1;
%         if f > c
%             break;
%         end
%     else 
%         f = 0;
%     end
%     Rx = [];
%     Ry = [];
%     for i = 1 : upxsize(1)
%         p = abs(upu(1)*upskedgex(i) + upu(2)*upskedgex(i+upxsize(1)) + upu(3) - upskedgey(i));
%         if p < updelta
%             Rx = [Rx; upskedgex(i), upskedgex(i+upxsize(1)), upskedgex(i+upxsize(1)*2)];
%             Ry = [Ry; upskedgey(i)];
%         end
%     end
%     uprx = Rx;
%     upry = Ry;
%     upr = uprxsize(1);
% end
% 
% upr = size(uprx);
% uplip = zeros(bbsize1, bbsize2, 1);
% for i = 1:upr(1)
%     uplip(upry(i), uprx(i+upr(1))) = 1;
% end
% uplip = double(uplip);
% % % figure, imshow(uplip);
% % % subplot(2, 3, 4), imshow(uplip);
% % figure, subplot(1, 3, 1), imshow(uplip);
% 
% f = 0;
% lowdelta = 1;
% lowr = 0;
% lowxsize = size(lowskedgex);
% lowrx = lowskedgex;
% lowry = lowskedgey;
% for t = 1 : 10
%     lowrxsize = size(lowrx);
%     lowx = pinv(lowrx);
%     lowu = lowx*lowry;
%     if lowrxsize <= lowr
%         f = f + 1;
%         if f > c
%             break;
%         end
%     else 
%         f = 0;
%     end
%     Rx = [];
%     Ry = [];
%     for i = 1 : lowxsize(1)
%         p = abs(lowu(1)*lowskedgex(i) + lowu(2)*lowskedgex(i+lowxsize(1)) + lowu(3) - lowskedgey(i));
%         if p < lowdelta
%             Rx = [Rx; lowskedgex(i), lowskedgex(i+lowxsize(1)), lowskedgex(i+lowxsize(1)*2)];
%             Ry = [Ry; lowskedgey(i)];
%         end
%     end
%     lowrx = Rx;
%     lowry = Ry;
%     lowr = lowrxsize(1);
% end
% lowr = size(lowrx);
% lowlip = zeros(bbsize1, bbsize2, 1);
% for i = 1:lowr(1)
%     lowlip(lowry(i), lowrx(i+lowr(1))) = 1;
% end
% lowlip = double(lowlip);
% % % figure, imshow(lowlip);
% % % subplot(2, 3, 5), imshow(lowlip);
% % subplot(1, 3, 2), imshow(lowlip);
% 
% edgepointsxy = [];
% for i = 1 : upr(1) + lowr(1)
%     if i <= upr(1)
%         edgepointsxy = [edgepointsxy; uprx(i+upr(1)), upry(i)];
%     elseif i >upr(1)
%         edgepointsxy = [edgepointsxy; lowrx(i-upr(1)+lowr(1)), lowry(i-upr(1))];
%     end
% end
% sizeofp = size(edgepointsxy);
% tr = zeros(bbsize1, bbsize2, 1);
% for i = 1:sizeofp(1)
%     tr(edgepointsxy(i+sizeofp(1)),edgepointsxy(i)) = 1;
% end
% tr = logical(tr);
% figure, imshow(tr);
% % subplot(2, 3, 6), imshow(tr);
% 
% % 找出BoundingBox的區域
% topwp = 0; t = 0;
% underwp = 0; u = 0;
% leftwp = 0; l = 0;
% rightwp = 0; r = 0;
% % find the top point
% for i = 1 : bbsize1   %mouthheight+1
%     for j = 1 : bbsize2   %mouthwidth+1
%         if tr(i,j) ~=0 && t == 0
%             topwp = i;
%             t = 1;
%         elseif t == 1
%             break;
%         end
%     end
% end
% % find the under point
% for i = bbsize1 : -1 : 1   %mouthheight+1 : -1 : 1
%     for j = bbsize2 : -1 : 1   %mouthwidth+1 : -1 : 1
%         if tr(i,j) ~=0 && u == 0
%             underwp = i;
%             u = 1;
%         elseif u == 1
%             break;
%         end
%     end
% end
% % find the left point
% for j = 1 : bbsize2   %mouthwidth+1
%     for i = 1 : bbsize1   %mouthheight+1
%         if tr(i,j) ~=0 && l == 0
%             leftwp = j;
%             l = 1;
%         elseif l == 1
%             break;
%         end
%     end
% end
% % find the right point
% for j =bbsize2 : -1 : 1   %mouthwidth+1 : -1 : 1
%     for i = bbsize1 : -1 : 1   %mouthheight+1 : -1 : 1
%         if tr(i,j) ~=0 && r == 0
%             rightwp = j;
%             r = 1;
%         elseif r == 1
%             break;
%         end
%     end
% end
% 
% wpwidth = rightwp - leftwp;
% wpheight = underwp - topwp;
% 
% rboundinglip = imcrop(lipROI, [leftwp topwp wpwidth wpheight]);
% 
% % % figure, subplot(1, 2, 1), imshow(boundinglip);
% % subplot(1, 3, 3), imshow(rboundinglip);
% % figure, imshow(tr);
% % figure, imshow(cannyedgeofROI);
% figure, imshow(rboundinglip);
% % imwrite(rboundinglip,resultdir);
% % 
% % graybb = rgb2gray(rboundinglip);
% % imwrite(graybb,graydir);
% % grayimg = 'newdatabase/data/1/ok/originalgray/user261001.jpg';
% % grayimg = imread(grayimg);
% 
% % xgray = double(grayimg);
% % mx = mean(xgray(:));
% % sx = std(xgray(:));
% % ygray = double(graybb);
% % my = mean(ygray(:));    % 算平均值
% % sy = std(ygray(:));     % 算標準差
% % 
% % normaly = (ygray-my)/sy*sx+mx;
% % normaly = uint8(normaly);
% % imwrite(normaly, nordir);
% % imwrite(normaly, lumdir);
% toc