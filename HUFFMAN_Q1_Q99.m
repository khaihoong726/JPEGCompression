clear all;
clc;

load M.mat;
load Lena_cropped.mat;

IMG_A1 = Lena10;

load Lena.mat;
I_OG = Lena;

IMG_full = Lena;

%% IN THIS CODE, WE WILL DO HUFFMAN ENCODING FOR QUANTIZATION MATRIX, Q1 TO Q99

compression_arr = [];
entropy_arr = [];
Eff_arr = [];
for compression = 0:10:100
    Q_matrix = quanmatrix(compression,M);
    fun = @(block_struct) round(dct2(block_struct.data)./Q_matrix);
    HUFF_Quant = blockproc(IMG_full,[8 8],fun, "PadPartialBlocks", true);
    symbols = unique(HUFF_Quant);
    counts = hist(HUFF_Quant(:), symbols);
    p = counts./sum(counts);
    [dict, avglen] = huffmandict(symbols,p);
    % Calculate the entropy
    p_arr = [];
    for i=1:length(p)
        p_arr(i) = p(i)*log2(1/p(i));  
    end
    entropy = sum(p_arr);
    % Calculate the efficiency of Huffman Encoding
    Efficiency = entropy/avglen*100;
    % APPEND THE ENTROPY OF ALL COMPRESSION LEVEL (1 TO 99) IN AN ARRAY
    compression_arr(end + 1) = compression;
    entropy_arr(end + 1) = entropy;
    % APPEND THE HUFFMAN CODE EFFICIENCY OF ALL COMPRESSION LEVEL (1 TO 99) IN AN ARRAY
    Eff_arr(end + 1) = Efficiency;
end

% PLOTTING ENTROPY, HUFFMAN EFFICIENCY AGAINST COMPRESSION LEVEL
figure(1)
subplot(2,1,1);
plot(compression_arr,entropy_arr),title("Entropy vs Q")
subplot(2,1,2);
plot(compression_arr,Eff_arr), title("Efficiency vs Q")

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