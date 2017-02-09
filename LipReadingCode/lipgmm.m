function [lipem, skinem] = lipgmm(lipcolordir,skincolordir)

inputlip_dir = dir(fullfile(lipcolordir,'*.jpg'));
[x, y] = size(inputlip_dir);
lipcolor = [];
for i = 1 : x
    lip = imread(fullfile(lipcolordir, inputlip_dir(i).name));
    lip = rgb2ycbcr(lip);
    szlip = size(lip);
    shlip = reshape(lip, szlip(1)*szlip(2),3);
    lipcolor = [lipcolor; shlip];
end
lipcolor = double(lipcolor);
lipem = gmdistribution.fit(lipcolor,1);

inputskin_dir = dir(fullfile(skincolordir,'*.jpg'));
[x, y] = size(inputskin_dir);
skincolor = [];
for i = 1 : x
    skin = imread(fullfile(skincolordir, inputskin_dir(i).name));
    skin = rgb2ycbcr(skin);
    szskin = size(skin);
    shskin = reshape(skin, szskin(1)*szskin(2),3);
    skincolor = [skincolor; shskin];
end
skincolor = double(skincolor);
skinem = gmdistribution.fit(skincolor,1);
end