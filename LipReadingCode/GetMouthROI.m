clear all;
tic;
% �إ߽���ήB��˥��ҫ�
lipcolordir = 'LipcolorSample';
skincolordir = 'SkincolorSample';
[firstlipem, firstskinem] = lipgmm(lipcolordir, skincolordir);

fid = fopen('fail.txt', 'wt');  % �إߤ@��.txt��(�˵L�k���L�B�ϰ쪺�ɮצW��)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};
for d = 1:1  % Ū�����O��Ƨ�
    cmddir = cmd{d};
for user = 1 : 90   % Ū���ϥΪ̸�Ƨ�
    userdir = num2str(user);
    datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\'];  % Ū���W�٬�'1-90'����Ƨ�
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI']);
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','lipROI']);
    input_dir = dir(fullfile(datadir, '*.jpg'));  % �N���w��Ƨ��������w���ɦW���Ҧ��ɮ׸�ơA��Jinput_dir�}�C
    [x, y ] = size(input_dir);
    for n = 1 : x
        [user,n]
        img = imread(fullfile(datadir, input_dir(n).name)); % Ū�����w��Ƨ��������w�ɮ�
        faceDetector =vision.CascadeObjectDetector;  % �H�y���Ѿ�
        fboxes = step(faceDetector, img);  % ��X�H�y
        %  fboxes[a,b,c,d]: a->���W����x�y�� b->���W����y�y�� c->�y���ϰ쪺�e d->�y���ϰ쪺��
        faceROI = insertObjectAnnotation(img, 'rectangle', fboxes, 'Face');
        [fboxesx, fboxesy] = size(fboxes);
        farea = 0; % �y�����n�Ѽ�
        maxarea = 0; % �̤j�y�����n�Ѽ�
        if fboxesx == 0 % �䤣��H�y�N�����ɮצW��
            fprintf(fid,'%s\n',input_dir(n).name);
            continue; % ���U�@�i�Ϥ�
        elseif fboxesx == 1 % ��X�y���ϰ쪺�I��e��
            left = fboxes(1);  % ���W���� x �y��
            top = fboxes(2) - 10;   % ���W���� y �y�� -10:�Ȥ����U�L�B����
            width = fboxes(3);  % �ϰ쪺�e
            height = fboxes(4) + 30;  % �ϰ쪺�� +30:�Ȥ����U�L�B����
        else % ��줣�u�@�i�y�ɡA���y�����n�̤j��
            for i = fboxesx*2+1:fboxesx*3 % �q�Ĥ@�Ӫ��e�����}�l
                farea = fboxes(i)*fboxes(i + fboxesx); % ���n=�e*��
                if farea > maxarea % �ثe�y���ϰ�j��̤j�y���ϰ�A�N���N�íק�y�лP�̤j���n��
                    maxarea = farea;
                    % ���y���ϰ쪺�y���I��e��
                    left = fboxes(i - 2*fboxesx);
                    top = fboxes(i - fboxesx) - 10;
                    width = fboxes(i);
                    height = fboxes(i + fboxesx) + 30;
                end
            end
        end
        faceROI = imcrop(img, [left top width height]);  % �N�H�y�ϰ�ŤU
        imwrite(faceROI,[datadir 'faceROI\' input_dir(n).name]);  % �N�y���ϰ�s���ɮ�
        
        mouthDetector = vision.CascadeObjectDetector('Mouth');  % �L�B���Ѿ�
        mboxes = step(mouthDetector, faceROI);  % ��X�L�B
        lipROI = insertObjectAnnotation(faceROI, 'rectangle', mboxes, 'Mouth');
        [mboxesx, mboxesy] = size(mboxes);
        mouthleft = 0;
        if mboxesx ==0 % �䤣��L�ڴN�����ɮצW��
            fprintf(fid,'%s\n',input_dir(n).name);
            continue; % ���U�@�i��
        end
        maxarea =2000;
        maxtop = 0; % ��̤j��Y
        re = 0;   % �P�_�O�_����ϰ쥢�Ѫ��Ѽ�
        if mboxesx == 1
            mouthleft = mboxes(1);
            mouthtop = mboxes(2) -10;
            mouthwidth = mboxes(3);
            mouthheight = mboxes(4) +15;
            re = 1; % �N�������\
        else
            for i = mboxesx+1:mboxesx*2  % ��Y��찵���P�_�з�
                larea = mboxes(i + mboxesx)*mboxes(i + mboxesx*2);  % ��L�ڰϰ쭱�n
                mouthright = mboxes(i - mboxesx) + mboxes(i + mboxesx);  % ��X�̤j��X
                %         if mboxes(i) > 2*height/3 && mboxes(i - mboxesx) < width/2 && mouthright > width/2 && larea > maxarea
                %     �L�B�ϰ�P�_���� : �̤U�����I�n>�y����6.3/10�A���k���䪺�I�n�b�y�����u�����k����A�i�઺�L�B�ϰ줤���n�̤j��
                if mboxes(i) > 6.3*height/10 && mboxes(i - mboxesx) < width/2 && mouthright > width/2 && larea > maxarea && mboxes(i) > maxtop
                    maxarea = larea;
                    maxtop = mboxes(i);
                    mouthleft = mboxes(i - mboxesx);
                    mouthtop = mboxes(i) - 10;
                    mouthwidth = mboxes(i + mboxesx);
                    mouthheight = mboxes(i + mboxesx*2) +15;
                    if mouthtop + mouthheight > height % �Y�O�L�B���̧C�I�W�X�y���A�h�H�y���̧C�I���L�B�̧C�I
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
        %  �ŤU�̫᪺�L�ڰϰ�æs��
        mouthroi = [mouthleft mouthtop mouthwidth mouthheight];
        lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);   % �ŤU�̫᪺�L�ڰϰ�
        imwrite(lipROI,[datadir 'mouthROI\' input_dir(n).name]);  % �N�L�ڰϰ�s���ɮ�
        %figure,imshow(lipROI);
        lipROIsize = size(lipROI);
        modle = zeros(lipROIsize(1), lipROIsize(2), 3, 'uint8');
        % �N�L�B�v���নYCbCr��m
        lipcolor = rgb2ycbcr(lipROI);
        slipcolor = size(lipcolor);
        % MxN���v�����s�ƦC��M*Nx1���}�C
        rlipcolor = reshape(lipcolor, slipcolor(1)*slipcolor(2),3);
        rlipcolor = double(rlipcolor);
        % �p��v�����C�ӹ����ݩ󽧦�ήB�⪺���v
        flipe = pdf(firstlipem,rlipcolor);
        fskine = pdf(firstskinem,rlipcolor);
        
        % �Q��YCbCr��X�L�B
        for p=1:slipcolor(1)*slipcolor(2)
            if (fskine(p)>flipe(p))
                rlipcolor(p,:) = [0.0627*255 0.5020*255 0.5020*255];
            end
        end
        rlipcolor = reshape(rlipcolor, slipcolor(1), slipcolor(2), 3);
        rlipcolor = uint8(rlipcolor);
        
        I = ycbcr2rgb(rlipcolor);  % ��^RGB�v��
        
        %figure, imshow(I);
        
        g = rgb2gray(I);  % RGB��Ƕ���
        bw = g > 0;  % �N�Ƕ��ϤG�Ȥ�
        
        % �G�ȹϰ��I�k
        se = strel('disk',4);
        erobw = imerode(bw,se);
        
        bcc = erobw;
        cc = bwconncomp(bcc);
        numpixels = cellfun(@length,cc.PixelIdxList);  % ��X�v����������(�զ�ϰ�)�����G
        [numsort, idxsort] = sort(numpixels);
        [biggest, bidx] = max(numpixels);
        second = 0; sid = 0;
        len = length(numsort);
        %figure, imshow(bcc);
        if len >= 2
            second = numsort(len-1);
            sid = idxsort(len-1);
        end
        %             �O�d�̤j�����զ�ϰ�
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
        
        % �Q�ΫI�k�᪺���G��X�s���B��ҫ��P�D�B��ҫ��إ�
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
        % �j���B�ⰻ��
        for i = 1 : 10
            testlip = lipROI;
            stestlip = size(testlip);
            rtestlip = reshape(testlip, stestlip(1)*stestlip(2),3);
            rtestlip = double(rtestlip);
            
            % pdf���v�ҫ� - �B��
            lipe = pdf(lipem,rtestlip);
            % pdf���v�ҫ� - ����
            skine = pdf(skinem,rtestlip);
            % ���e��theta
            theta = [5e-6];
            
            % �B��P�w
            for p=1:stestlip(1)*stestlip(2)
                if (skine(p)>lipe(p))
                    rtestlip(p,:) = [0 0 0];
                end
            end
            
            result  = reshape(rtestlip, stestlip(1), stestlip(2), 3);
            result = uint8(result);
            
            testliparray = []; % �s���B��˥��}�C
            testskinarray = []; % �s������˥��}�C
            
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
            lipem = gmdistribution.fit(testliparray,2); % �s���B��ҫ�
            
            testskinarray = double(testskinarray);
            skinem = gmdistribution.fit(testskinarray,2); % �s������ҫ�
            
            sumlip = size(testliparray);
            sumlipcolorarr = [sumlipcolorarr; sumlip(1)];
            
            if i >=3
                if sumlipcolorarr(i) - sumlipcolorarr(i-1) < 0.78*(sumlipcolorarr(i-1) - sumlipcolorarr(i-2))
                    break;
                end
            end
        end
        
        %figure, imshow(result);
        
        % �j�馡���G�I�k
        g = rgb2gray(result);  % RGB��Ƕ���
        bw = g > 0;  % �N�Ƕ��ϤG�Ȥ�
        
        % �G�ȹϰ�closing
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
        
        % ��XBoundingBox���ϰ�
        topwp = 0; t = 0;
        underwp = 0; u = 0;
        leftwp = 0; l = 0;
        rightwp = 0; r = 0;
        % ��X�̤W�����I
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
        % ��X�̤U�����I
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
        % ��X�̥��䪺�I
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
        % ��X�̥k�䪺�I
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
        % ��boundingbox�����B���G�����u���X�Ay = ax^2+bx+c�A���P�_�����O�_���B��A�Y�ݩ�B��h�d�U
        for i = 1 : bbsize1
            for j = 1 : bbsize2
                points = [lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
                points = double(points);
                % �p��BoundingBox���������ݩ�B�⪺���v
                flipe = pdf(lipem, points);
                fskine = pdf(skinem, points);
                if bwtrace(i,j) == 1 && i <= bbsize1/2 % �W�B
                    if fskine < flipe
                        %                 upskedge = [upskedge; lipROI(i,j,1), lipROI(i,j,2), lipROI(i,j,3)];
                        upskedgey = [upskedgey; i];                           % y = [y1;y2;...;yn]
                        upskedgex = [upskedgex; j^2, j, 1];              % x = [x1^2, x1, 1; x2^2, x2, 1; ...; xn^2, xn, 1]
                    else
                        bwtrace(i,j) = 0;
                    end
                elseif bwtrace(i,j) == 1 && i > bbsize1/2 % �U�B
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
        
        % �P�_�ݩ�B�⪺�I�O�_�����u�W���I-�W�B
        f = 0;
        upr = 0;
        c = 3;
        updelta = 3;
        upxsize = size(upskedgex);    % Eskin���ƶq
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
        
        % �P�_�ݩ�B�⪺�I�O�_�����u�W���I-�U�B
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
        
        % ��XBoundingBox���ϰ�
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
        % �ťX�̫᪺�L�B�ϰ�
        rboundinglip = imcrop(lipROI, [leftwp topwp wpwidth wpheight]);
        imwrite(rboundinglip,[datadir 'lipROI\' input_dir(n).name]);
        %figure, imshow(rboundinglip);
%         %             imwrite(rboundinglip,[lipdir '\result\' inputlip_dir(n).name]);
%         imwrite(rboundinglip,[datadir '\result\' input_dir(n).name]);
        % �ন�Ƕ��v��
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
        %             my = mean(ygray(:));    % �⥭����
        %             sy = std(ygray(:));     % ��зǮt
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