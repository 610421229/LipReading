clear all;
tic;

fid = fopen('fail.txt', 'wt');  % �إߤ@��.txt��(�˵L�k���L�B�ϰ쪺�ɮצW��)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};

for d = 2:6  % Ū�����O��Ƨ�
    
    cmddir = cmd{d};
    
    for user = 1 : 90   % Ū���ϥΪ̸�Ƨ�
        
        userdir = num2str(user);
        datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\'];  % Ū���W�٬�'1-90'����Ƨ�
        mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','faceROI']);
        mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
        input_dir = dir(fullfile(datadir, '*.jpg'));  % �N���w��Ƨ��������w���ɦW���Ҧ��ɮ׸�ơA��Jinput_dir�}�C
        [x, y ] = size(input_dir);
        
        for n = 1 : x
                
            [user,n]
            img = imread(fullfile(datadir, input_dir(n).name)); % Ū�����w��Ƨ��������w�ɮ�
            img = imresize(img,[360 480]);
            faceDetector =vision.CascadeObjectDetector;  % �H�y���Ѿ�
            fboxes = step(faceDetector, img);  % ��X�H�y
            %  fboxes[a,b,c,d]: a->���W����x�y�� b->���W����y�y�� c->�y���ϰ쪺�e d->�y���ϰ쪺��
            faceROI = insertObjectAnnotation(img, 'rectangle', fboxes, 'Face');
            [fboxesx, fboxesy] = size(fboxes);
            farea = 0; % �y�����n�Ѽ�
            maxarea = 0; % �̤j�y�����n�Ѽ�
            
                if fboxesx == 0 % �䤣��H�y�N�����ɮצW��
                    %fprintf(fid,'%s\n',input_dir(n).name);
                    %continue; % ���U�@�i�Ϥ�
                
                elseif fboxesx == 1 % ��X�y���ϰ쪺�I��e��
                    
                    left = fboxes(1);  % ���W���� x �y��
                    top = fboxes(2) - 10;   % ���W���� y �y�� -10:�Ȥ����U�L�B����
                    width = fboxes(3);  % �ϰ쪺�e
                    height = fboxes(4) + 30;  % �ϰ쪺�� +30:�Ȥ����U�L�B����
                    
                else % ��줣�u�@�i�y�ɡA���y�����n�̤j��
                    
                    for i = fboxesx*2+1 : fboxesx*3 % �q�Ĥ@�Ӫ��e�����}�l
                        
                            farea = fboxes(i)*fboxes(i + fboxesx); % ���n=�e*��
                        
                            if farea > maxarea % �ثe�y���ϰ�j��̤j�y���ϰ�A�N���N�íק�y�лP�̤j���n��
                            
                            maxarea = farea;
                            % ���y���ϰ쪺�y���I��e��
                            left = fboxes(i - 2*fboxesx);
                            top = fboxes(i - fboxesx) - 10;
                            width = fboxes(i);
                            height = fboxes(i + fboxesx) ;
                            
                            end
                    end
                end
        faceROI = imcrop(img, [left top+height*2/3 width height/3]);% �N�H�y�ϰ�ŤU
        %faceROI = imcrop(faceROI, [left top+height*2/3 width height]);
        imwrite(faceROI,[datadir 'faceROI\' input_dir(n).name]);  % �N�y���ϰ�s���ɮ�
        
        mouthDetector = vision.CascadeObjectDetector('Mouth');  % �L�B���Ѿ�
        mboxes = step(mouthDetector, faceROI);  % ��X�L�B
        lipROI = insertObjectAnnotation(faceROI, 'rectangle', mboxes, 'Mouth');
        [mboxesx, mboxesy] = size(mboxes);
        mouthleft = 0;
        maxarea =2000;
        maxtop = 0; % ��̤j��Y
        re = 0;   % �P�_�O�_����ϰ쥢�Ѫ��Ѽ�
        
        if mboxesx ==0 % �䤣��L�ڴN�����ɮצW��
            lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);
            %fprintf(fid,'%s\n',input_dir(n).name);
            %continue; % ���U�@�i��
            
        elseif mboxesx == 1
            mouthleft = mboxes(1);
            mouthtop = mboxes(2) -10;
            mouthwidth = mboxes(3);
            mouthheight = mboxes(4) +15;
            re = 1; % �N�������\
        else

            for i = mboxesx+1:mboxesx*2  % ��Y��찵���P�_�з�
                larea = mboxes(i + mboxesx)*mboxes(i + mboxesx*2);  % ��L�ڰϰ쭱�n
            
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
        lipROI = imcrop(faceROI, [mouthleft mouthtop mouthwidth mouthheight]);   % �ŤU�̫᪺�L�ڰϰ�
        imwrite(lipROI,[datadir 'mouthROI\' input_dir(n).name]);  % �N�L�ڰϰ�s���ɮ�

        end
    end
end
fclose(fid);
toc;