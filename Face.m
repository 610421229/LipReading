clear all;
tic;
%% �إ߽���ήB��˥��ҫ�
% lipcolordir = 'LipcolorSample';
% skincolordir = 'SkincolorSample';
% [firstlipem, firstskinem] = lipgmm(lipcolordir, skincolordir);
%%
fid = fopen('fail.txt', 'wt');  % �إߤ@��.txt��(�˵L�k���L�B�ϰ쪺�ɮצW��)

cmd = {'Drink' 'Eat' 'Spa' 'Walk' 'Shower' 'Toilet'};
for d = 1:1  % Ū�����O��Ƨ�
    cmddir = cmd{d};
for user = 1 : 1   % Ū���ϥΪ̸�Ƨ�
    userdir = num2str(user);
    datadir = ['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\'];  % Ū���W�٬�'1-90'����Ƨ�
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','mouthROI']);
    mkdir(['MatlabDataBase\six_cmd\images\' cmddir '\' userdir '\','lipROI']);
    input_dir = dir(fullfile(datadir, '*.jpg'));  % �N���w��Ƨ��������w���ɦW���Ҧ��ɮ׸�ơA��Jinput_dir�}�C
    [x, y ] = size(input_dir);
     for n = 28 : 29
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