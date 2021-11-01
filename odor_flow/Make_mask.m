%Make_mask

%% load image
M = imread('Z:\Kevin\20190817\Data20190817_165031\Frame_000000.jpg');
sx = size(M,1);
sy = size(M,2);
imagesc(M)

%% filling space
center = [sx/2, sy-1250];
diam = 1300;
Mask = zeros(size(M));
for ii = 1:sx
    for jj = 1:sy
        if sqrt(sum(([ii,jj]-center).^2)) < diam
            Mask(ii,jj) = 1;
        end
    end
end
figure;
imagesc(Mask)

%% binerize  (binarize to [0,255] image)
pos1 = find(Mask==0);
Mask(pos1) = 255;
pos2 = find(Mask==1);
Mask(pos2) = 0;
figure;
imagesc(Mask)

%% Save (save as tif file)
f_name = 'Flow_Mask.tif';
imwrite(Mask,f_name);


