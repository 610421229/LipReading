clc
clear

cmd = 'Walk';
inputdir = ['MatlabDataBase/six_cmd/videos/' cmd '/']; %�v�����|��Ƨ�
outputd = ['MatlabDataBase/six_cmd/images/' cmd '/']; %�Ϥ���X��Ƨ�
input_dir = dir(fullfile(inputdir, '*.mp4')); %Ū�����b��Ƨ��U�Ҧ�mp4�ɪ��ԲӸ�T
[x, y] = size(input_dir); %x���ɮ׭Ӽ�
for n = 1 : 1
n
    videonum = num2str(n); % ��e�j�����r��
    outputdir = [outputd videonum '/'];
    mkdir(outputdir); % �s�W�@�ӿ�X��Ƨ�
    obj = VideoReader(fullfile(inputdir, input_dir(n).name)); %�Ϩ���e�v���ԲӸ�T
    objn = obj.NumberOfFrames; % Ū����e�v��Frame��
    objw = obj.Width; %�v���e
    objh = obj.Height; %�v����

    mov(1:objn) = struct('cdata',zeros(objw,objh,3,'uint8'),'colormap',[]); %�غccdata�s720*960*3 Unit8�Ϥ� ,colormap�s�ůx�}
    a = size(input_dir(n).name); %���o�v���ɮצW�١A�Ҧp : Walk(1).mp4
    d = a(2)-4; %�ťh�ɦW4�̫�Ӧr(.mp4)
    f = input_dir(n).name(1:d); %�s���ѤU�ɦW(Walk(1))
    for k=1:objn
        mov(k).cdata = read(obj,k);
        output = mov(k).cdata;  %�s����e�Ϥ�
        tempStr = strcat(num2str(k),'.jpg'); % �s���ɦW�A�Ҧp:1.jpg

        if k < 10
            tempStr = [outputdir f '00' tempStr]; % �N�ɦW�p��10����00�A�Ҧp : 001.jpg
        else
            tempStr = [outputdir f  '0' tempStr]; %�N�ɦW�j��10����0�A�Ҧp : 010.jpg
        end
        imwrite(output,tempStr); %�N��e�ɮ׼g�J���w���|
    end
end