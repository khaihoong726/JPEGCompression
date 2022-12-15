clear all;
clc;

load M.mat;
load Lena_cropped.mat;

IMG_A1 = Lena10;

load Lena.mat;
I_OG = Lena;

M %Printing Matrix of Pixel Values of the ROI

% IN THIS SECTION, WE WILL DO HUFFMAN CODING AFTER ZIGZAG RLE APPROACH 
% TOWARDS THE FULL 8 X 8 CROPPED LENA IMAGE. PREVIOUSLY, WE DO HUFFMAN CODING
% ON WHOLE 512 X 512 DIRECTLY. 

%VALUES SHIFTED TO CENTRE AROUND 0

for i=1:length(M)
    for j=1:length(M)
        IMG_A2(i,j) = IMG_A1(i,j)-128;
    end
end
IMG_A2
compression = 90; % Input desired Quantization Matrix, Q
Qn = quanmatrix(compression,M)
Q50 = M;
dctIMG_A = dct2(IMG_A2);
QUANT = round(dctIMG_A./Qn)

% ZIGZAG RLE PATTERN
ZZ = zigzag(QUANT)
IZZ = izigzag(ZZ,8,8)

% HUFFMAN + ZIGZAG
Tcount = size(ZZ,1)*size(ZZ,2);
symbol = unique(ZZ); % ENCODING THE VECTORS FROM ZIGZAG RLE
counts = histc(ZZ(:), symbol);
P = counts./Tcount;
% DECLARING DICTIONARY FOR AVG CODEWORD FOR EACH SYMBOL FOR HUFFMAN CODING
[dict, avglen] = huffmandict(symbol, P);
% RESHAPE THE MATRIX 
new_shape = reshape(ZZ.',1,[]);
% ENCODE USING HUFFMAN
HuffmanCode = huffmanenco(new_shape,dict);

H = [];
for i=1:length(P)
    for j=1:length(symbol)
        H(i) = P(i)*log2(1/P(i));
    end
end
H_Huffman = sum(H);
% CALCULATING EFFICIENCY OF HUFFMAN CODES
avglen
H_Huffman
Efficiency = H_Huffman/avglen*100

% DECODE USING HUFFMAN
HuffmanDecode = huffmandeco(HuffmanCode, dict);
Huffman_IZZ = izigzag(HuffmanDecode,8,8)

[r,c] = size(I_OG);
centering(1:r,1:c) = 128;
Lena_c = I_OG - centering;

% HERE, WE WANT TO COMPARE IF THERE IS ANY DIFFERENCES IN DOING ZIG-ZAG
% PATTERN RLE TO THE IMAGE

% W/O ZIG ZAG RLE
fun1 = @(matrix) round((dct2(matrix))./Qn);
fun2 = @(matrix) idct2(matrix.*Qn);
dct_C = blkproc(Lena_c,[8 8],fun1);
idct_c = blkproc(dct_C,[8 8],fun2);
idct_c = idct_c + centering;

% WITH ZIG ZAG RLE
zigzagfun = @(matrix) zigzag(matrix);
izigzagfun = @(matrix) izigzag(matrix,8,8);
Zigzag_block = blkproc(dct_C,[8 8],zigzagfun);

% HUFFMAN AFTER ZIGZAG
ZZCount = size(Zigzag_block,1)*size(Zigzag_block,2);
ZZSym = unique(Zigzag_block);
ZZCounts = histc(Zigzag_block(:), ZZSym);
P_ZZ = ZZCounts./ZZCount;
[ZZ_dict, ZZ_avglen] = huffmandict(ZZSym, P_ZZ);
ZZ_reshape = reshape(Zigzag_block,1,[]);
ZZ_HUFF = huffmanenco(ZZ_reshape,ZZ_dict);

% HUFFMAN EFFICIENCY OF ZIGZAG
H_1 = [];
for i=1:length(P_ZZ)
    for j=1:length(ZZSym)
        H_1(i) = P_ZZ(i)*log2(1/P_ZZ(i));
    end
end
H_ZZ = sum(H_1);
Eff_ZZ = H_ZZ/ZZ_avglen*100

% HUFFMAN DECODE
ZZ_HUFFDE = huffmandeco(ZZ_HUFF,ZZ_dict);
ZZ_reshape_DE = reshape(ZZ_HUFFDE,size(Zigzag_block));
Izigzag_block = blkproc(ZZ_reshape_DE,[1 64], izigzagfun);
Zigzag_IDCT = blkproc(Izigzag_block,[8 8],fun2);
Zigzag_IDCT = Zigzag_IDCT + centering;

figure(1);
subplot(1,2,1), imshow(Zigzag_IDCT,[]), title("Q" + compression +  " with Zig Zag", FontSize=18);
subplot(1,2,2), imshow(idct_c,[]), title("Q" + compression +  " W/O Zig Zag",FontSize=18);

immse(Zigzag_IDCT,idct_c)

%% FUNCTION TO CALCULATE DIFFERENT COMPRESSION LEVEL Q MATRIX
function QM = quanmatrix(F,M)
    if F >= 50
        S = 200 - 2*F;
    else 
        S = 5000/F;
    end 
    QM = floor((50+S.*M)/100);
    % PREVENT DIVISION OF 0
    QM(QM==0) = 1;
    QM(QM>=255) = 255;
end
