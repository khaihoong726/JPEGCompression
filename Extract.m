load lena.mat
img = Lena;
figure(1)
imshow(img, [])
targetSize = [8 8];
%imcrop
% Position of Lena's Right Eye, cropped with the size of 8x8 
% for [x, y, n-1, n-1]
Lena10 = imcrop(img,[325.5 262.5 7 7]) 

figure(2)
subplot(2,1,1);
imshow(Lena10,[]);
title('Lena Eye')