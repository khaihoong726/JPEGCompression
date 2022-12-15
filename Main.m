%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         PART A            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
clc;

% LOADING LENA 
load Lena.mat;
I_OG = Lena;

% SHOWING ORIGINAL LENA IMAGE
figure(1);
imshow(I_OG,[])
title('Original 512x512 Lena Image'); % Default is grey scale

%%
% Show ROI of 8x8 from Lena Image
load M.mat;
load Lena_cropped.mat;

IMG_A1 = Lena10;
figure(2);
subplot(1,2,2);
imshow(IMG_A1,[])
title('8x8 pixels (W/O Quantization)'); % Default is grey scale

M %Printing Matrix of Pixel Values of the ROI

%VALUES SHIFTED TO CENTRE AROUND 0

for i=1:length(M)
    for j=1:length(M)
        IMG_A2(i,j) = IMG_A1(i,j)-128;
    end
end
IMG_A2

%%
% PERFORM DCT
dctIMG_A = dct2(IMG_A2)

%%
% PERFORM IDCT
idctIMG_A = idct2(dctIMG_A)

%%
% VALUES SHFITED BACK TO ORIGINAL OFFSET
for i=1:length(M)
    for j=1:length(M)
        idctIMG_A2(i,j) = idctIMG_A(i,j) + 128;
    end 
end 
Result_A = idctIMG_A2
subplot(1,2,1);
imshow(Result_A,[]);
title('Compressed 8x8 pixels (w/o Quantization')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         PART B            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load M.mat; % LOADING JPEG DEFAULT QUANTIZATION MATRIX Q50

Q50 = M;

% QUANTIZATION WITH Q50
% DIVIDE THE DCT IMAGE BY THE QUANTIZATION MATRIX AND ROUNDING THE RESULT
% TO THE NEAREST INTEGER
QUANT = round(dctIMG_A./Q50)

%%
% DEQUANTIZATION
% MULTIPLY BY THE QUANTIZATION MATRIX TO UNDO THE QUANTIZATION

DQUANT = round(QUANT.*Q50)

%%
% PERFORM INVERSE DCT
 idctIMG_B = idct2(DQUANT);
for i = 1:8
    for j= 1:8
        idctIMG_B2(i,j) = idctIMG_B(i,j) + 128;
    end
end

Result_B = idctIMG_B2

%% 
% COMPARISON OF RESTORED PIXEL VALUES WITH ORIGINAL VALUE
figure(3)
subplot(2,1,1);
imshow(IMG_A1,[]);
title('Original 8x8 Crop of Lena''s Right Eye','FontSize', 30);
subplot(2,1,2);
imshow(Result_B,[]);
title('Compressed 8x8 Pixels (Q50)', 'FontSize', 30)

imwrite(uint8(Result_B), 'Lena_Q50.jpeg');
imwrite(uint8(IMG_A1), 'Lena_Cropped_ori.jpeg');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         PART C            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% CENTRALISE FULL 512X512 PIXELS IMAGE

for i=1:length(I_OG)
    for j=1:length(I_OG)
        I_OGCENTERED(i,j) = I_OG(i,j) - 128;
    end
end

I_OGCENTERED

%%
% PERFORM DCT ON THE FULL 512 X 512 PIXELS IMAGE

dctFull = dct2(I_OGCENTERED);

dctFull

%% 
% PERFORM QUANTIZATION OF THE FULL IMAGE WITH Q50
% FOR FULL IMAGE, DIVIDE THE BLOCKS OF 8 X 8 PIXELS IMAGE BY Q50 MATRIX AND
% ROUNDING IT OFF TO NEAREST INTEGER

QUANTFULL = @(block_struct) round(block_struct.data./Q50);
FULL_IMG_QUANTFULL = blockproc(dctFull, [8 8], QUANTFULL);

%% 
% PERFORM DEQUANTIZATION OF THE FULL IMAGE WITH Q50
% FOR FULL IMAGE, MULTIPLY THE BLOCKS OF 8 X 8 PIXELS IMAGE BY Q50 MATRIX AND
% ROUNDING IT OFF TO NEAREST INTEGER

DQUANTFULL = @(block_struct) round(block_struct.data.*Q50);
FULL_IMG_DQUANTFULL = blockproc(FULL_IMG_QUANTFULL, [8 8], DQUANTFULL);

%%
% PEFORM INVERSE DCT
idctFull = idct2(FULL_IMG_DQUANTFULL);

%%
% RETURNING THE VALUES OF 512 X 512 PIXELS TO ORIGINAL
for i=1:length(I_OG)
    for j=1:length(I_OG)
        idctFull2(i,j) = idctFull(i,j) + 128;
    end
end

idctFull2;

figure(4)
subplot(1,2,1);
imshow(I_OG,[]);
title('Original 512x512 Pixels');
subplot(1,2,2);
imshow(idctFull2,[])
title('Compressed 512 x 512 Pixels at Q50')

%%
% DOING FOR 10 DIFFERENT QUALITY FACTORS
% For 10 different quality factors

[r,c] = size(I_OG);
centering(1:r,1:c) = 128;
Lena_c = I_OG - centering;
b= 1; mse = []; psnr = []; 
% IDCT_mat = {}; DCT_mat = {}; Q = {};

% FOR Q10 TO Q100 AT INCREMENT OF 10
for i= 10:10:100
    Q = quanmatrix(i,M);
    fun1 = @(matrix) round((dct2(matrix))./Q);
    fun2 = @(matrix) idct2(matrix.*Q);
    dct_C = blkproc(Lena_c,[8 8],fun1);
    idct_c = blkproc(dct_C,[8 8],fun2);
    idct_c = idct_c + centering;
    if rem((i/10),2) ~= 0 || i == 100
%         DCT_mat{b} = dct_C;
%         IDCT_mat{b} = idct_c; 
%         Q_mat{b} = Q;
        figure(5), subplot(2,3,b), imshow(idct_c,[]);
        if i == 100
            title("No Compression");
        else
            title("Q" + i);
        end 
        mse(b) = immse(idct_c, I_OG);
        psnr(b) = 10*log10((255^2)/mse(b));
        B = idct_c./255;
        imwrite(B,['Q ', num2str(i), '%.tif']);
        b = b + 1;
    end 

    % CHANGE i TO ALTER THE QUALITY FACTOR, Q TO BE ENCODED USING HUFFMAN,
    % i.e. Q10 ==> i == 10
    if i == 10 % QUALITY THAT YOU WANT TO SAVE
        % SAVE FOR THE USE OF HUFFMAN CODING LATER
        Q100 = Q;
        Q100_dct = dct_C;
        Q100_idct = idct_c;
    end
end 

% SHOWING DATA REDUCTION RATE
zero_Q = 0;
non_zero = 0;

% COUNTING THE NUMBER OF ZERO ELEMENTS AFTER QUANTIZATION
for R = 1:size(Q100_dct,1)
    for C = 1:size(Q100_dct,2)
        if Q100_dct(R,C) == 0
            zero_Q = zero_Q + 1;
        end
    end
end

% COUNTING THE NUMBER OF NON-ZERO ELEMENTS IN ORIGINAL LENA
for R = 1:size(Q100_dct,1)
    for C = 1:size(Q100_dct,2)
        if Lena(R,C) ~= 0
            non_zero = non_zero + 1;
        end
    end
end

non_zero
zero_Q

% NON-ZERO ELEMENTS IN THE QUANTIZED MATRIX
Non_zero = non_zero - zero_Q;
% FORMULA TO GET DATA REDUCTION RATE
Data_Reduction = (1-Non_zero/size(Q100_dct,1)/size(Q100_dct,2))*100

% PLOTTING OF 
Q10 = quanmatrix(10,M);
compression_level = {'Q10', 'Q30', 'Q50', 'Q70', 'Q90'}';
Table_C = table(compression_level, mse(1:5)', psnr(1:5)' , 'VariableNames',{'Compression Level', 'MSE', 'PSNR (in dB)'})

figure(6);
plot(10:20:90,mse(1:5));
grid on;
title("Mean Squared Error in Compression", FontSize=25);
xlabel("Quality Factor, Q", FontSize=18);
ylabel("MSE", FontSize=18)

func1 = @(matrix) round((dct2(matrix))./Q10);
func2 = @(matrix) idct2(matrix.*Q);
dct_C2 = blkproc(Lena_c,[8 8],func1); 
idct_C2 = blkproc(dct_C2,[8 8],fun2);

figure(7);
% DISPLAY THE FINAL IMAGE OF LENA WITH QUANTIZATION MATRIX SAVED IN LINE 201
imshow(idct_C2,[]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         PART D            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eff = [];
% 
% for i=1:size(DCT_mat,2)
%     [m,n] = size(DCT_mat{1,i});
%     Tcount = m*n;
%     symbol = unique(DCT_mat{1,i});
%     counts = histc(DCT_mat{1,i}(:), symbol);
%     P = counts./Tcount;
%     % DECLARING DICTIONARY FOR AVG CODEWORD FOR EACH SYMBOL FOR HUFFMAN CODING
%     [dict, avglen] = huffmandict(symbol, P);
%     % RESHAPE THE MATRIX 
%     new_shape = reshape(DCT_mat{1,i}.',1,[]);
%     % ENCODE USING HUFFMAN
%     HuffmanCode = huffmanenco(new_shape,dict);
%     % DECODE USING HUFFMAN
%     HuffmanDecode = huffmandeco(HuffmanCode, dict);
%     % RESHAPE DECODED SIGNAL TO 512 X 512
%     Q100_Decode = reshape(HuffmanDecode,512,512);
%     % PERFORM BLOCK IDCT TO THE DECODED RESHAPED SIGNAL
%     fun2 = @(matrix) idct2(matrix.*Q_mat{1,i});
%     Q100D_idct = blkproc(Q100_Decode,[8 8],fun2);
%     Q100D_idct = Q100D_idct + centering;
%     % CALCULATING THE ERROR BETWEEN HUFFMAN ENCODED AND DECODED SIGNAL
%     MSE_HUFF = immse(DCT_mat{1,i}, Q100D_idct.')
%     H = [];
%     for k=1:length(P)
%         for l=1:length(symbol)
%             H(k) = P(k)*log2(1/P(k));
%         end
%     end
%     H_Huffman = sum(H);
%     Eff(i) = H_Huffman/avglen*100;
% end

[m,n] = size(Q100_dct);
Tcount = m*n;
symbol = unique(Q100_dct);
counts = histc(Q100_dct(:), symbol);
P = counts./Tcount;
% DECLARING DICTIONARY FOR AVG CODEWORD FOR EACH SYMBOL FOR HUFFMAN CODING
[dict, avglen] = huffmandict(symbol, P);
% RESHAPE THE MATRIX 
new_shape = reshape(Q100_dct.',1,[]);
% ENCODE USING HUFFMAN
HuffmanCode = huffmanenco(new_shape,dict);
% DECODE USING HUFFMAN
HuffmanDecode = huffmandeco(HuffmanCode, dict);
% RESHAPE DECODED SIGNAL TO 512 X 512
Q100_Decode = reshape(HuffmanDecode,512,512);
% PERFORM BLOCK IDCT TO THE DECODED RESHAPED SIGNAL
fun2 = @(matrix) idct2(matrix.*Q100);
Q100D_idct = blkproc(Q100_Decode,[8 8],fun2);
Q100D_idct = Q100D_idct + centering;
% CALCULATING THE ERROR BETWEEN HUFFMAN ENCODED AND DECODED SIGNAL
MSE_HUFF = immse(Q100_idct, Q100D_idct.')

% CALCULATING ENTROPY OF HUFFMAN CODING
H = [];
for i=1:length(P)
    for j=1:length(symbol)
        H(i) = P(i)*log2(1/P(i));
    end
end
H_Huffman = sum(H);
% CALCULATING EFFICIENCY OF HUFFMAN CODES
Efficiency = H_Huffman/avglen*100

figure(8);
subplot(1,2,1), imshow(Q100_idct,[]), title('Huffman Encoded');
subplot(1,2,2), imshow(Q100D_idct.',[]), title('Huffman Decoded');

imwrite(uint8(Q100_idct), 'Huffman_Encoded.jpeg');
imwrite(uint8(Q100D_idct.'), 'Huffman_Decoded.jpeg');
imwrite(uint8(Q100_idct), 'Q90_Lena.jpeg');
imwrite(uint8(Lena), 'Original_Lena.jpeg');

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
