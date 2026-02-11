function score = DTD(imageRef, imageDis)
% ========================================================================
% DTD Index with automatic downsampling, Version 1.0 2026.2.11
% Copyright(c) 2026 Chenyang Shi
% All Rights Reserved.
%史晨阳*，吴俊杰，袁瀚成，吴路路. 双尺度画面细节信息引导的全参考图像质量客观评价, 
%光学精密工程，2026.
%SHI Chenyang*,Wu Junjie, Yuan Hancheng, et al.Dual-scale tableau detail information-guided 
%full-reference objective image quality assessment[J].Optics and Precision Engineering,
%2026.(in chinese)
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------


if isequal(imageRef, imageDis)
    score = 0;
    return;
end

%% 颜色空间转换（同时生成亮度+色度）
[Ref_L, Ref_H, Ref_M] = rgb2opponent_full(imageRef); % 修改函数
[Dist_L, Dist_H, Dist_M] = rgb2opponent_full(imageDis);


%% 动态下采样（所有通道同步处理）
[Ref_L, Ref_H, Ref_M] = mdsi_downsample(Ref_L, Ref_H, Ref_M); % 参数修改
[Dist_L, Dist_H, Dist_M] = mdsi_downsample(Dist_L, Dist_H, Dist_M);

dx = [1 0 -1; 1 0 -1; 1 0 -1]/3;
dy = dx';
C1 = 140; 
C3 = 550;

   
    
    [gR, gD, gF] = compute_gradient_maps(Ref_L, Dist_L,dx, dy);
    GS12 = (2*gR.*gD + C1) ./ (gR.^2 + gD.^2 + C1);
    GS13 = (2*gR.*gF + C1) ./ (gR.^2 + gF.^2 + C1);
    GS23 = (2*gD.*gF + C1) ./ (gD.^2 + gF.^2 + C1);
    gradSim = GS12 + GS23 - GS13;

    % 保持原颜色相似性计算
    colorSim = enhanced_color_similarity(Ref_H, Ref_M, Dist_H, Dist_M, C3);
    
    % 画面细节累积
    saliencyMap1 = visualSaliency(Ref_L);
    saliencyMap2 = visualSaliency(Dist_L);
     weightedDiff =saliencyMap1 + saliencyMap2;  
     saliencyWeight=weightedDiff;

    combinedSim = 0.7*std2(gradSim) + 0.25*std2(colorSim)+0.05*std2(saliencyWeight);

 score=combinedSim;
end

function [gR, gD, gF] = compute_gradient_maps(Ref_L,Dist_L,dx, dy)
% 输入:
%   Ref_L: 参考图像亮度通道矩阵
%   Dist_L: 失真图像亮度通道矩阵
%   dx: Prewitt水平梯度算子
%   dy: Prewitt垂直梯度算子
% 输出:
%   gR, gD, gF: 参考、失真、融合图像的梯度幅值图

% Step 1: 生成融合亮度图（MDSI论文式(3)）
F =  0.5*Ref_L+0.5*Dist_L ;

% Step 2: 计算参考图像梯度
IxR = imfilter(Ref_L, dx, 'same', 'replicate'); % 水平梯度
IyR = imfilter(Ref_L, dy, 'same', 'replicate'); % 垂直梯度
gR = sqrt(IxR.^2 + IyR.^2);                    % 梯度幅值（L2范数）

% Step 3: 计算失真图像梯度
IxD = imfilter(Dist_L, dx, 'same', 'replicate');
IyD = imfilter(Dist_L, dy, 'same', 'replicate');
gD = sqrt(IxD.^2 + IyD.^2);


% Step 4: 计算融合图像梯度
IxF = imfilter(F, dx, 'same', 'replicate');
IyF = imfilter(F, dy, 'same', 'replicate');
gF = sqrt(IxF.^2 + IyF.^2);
end

%% 改进的颜色空间转换函数
function [L, H, M] = rgb2opponent_full(img)
% 同时生成亮度L和色度H/M
R = double(img(:,:,1));
G = double(img(:,:,2));
B = double(img(:,:,3));
%%LHM颜色通道
L = 0.2989*R + 0.5870*G + 0.1140*B; % 标准亮度计算
H = 0.30*R + 0.04*G - 0.35*B;      % 色度H通道
M = 0.34*R - 0.60*G + 0.17*B;      % 色度M通道

%%YUV颜色通道
% L = 0.299*R + 0.587*G + 0.114*B; % 标准亮度计算
% H = -0.147*R - 0.289*G + 0.436*B;      % 色度H通道
% M = 0.615*R - 0.515*G - 0.100*B;      % 色度M通道

% %%YIQ颜色通道
% L = 0.299*R + 0.587*G + 0.114*B; % 标准亮度计算
% H =  0.596*R - 0.274*G - 0.322*B;      % 色度H通道
% M = 0.211*R - 0.523*G + 0.312*B;      % 色度M通道

%%LMN颜色通道
% L = 0.06*R + 0.63*G + 0.27*B; % 标准亮度计算
% H = 0.30*R + 0.04*G - 0.35*B;      % 色度H通道
% M = 0.34*R - 0.60*G + 0.17*B;      % 色度M通道

%%YCbCr颜色通道
% L = 0.2126*R + 0.7152*G + 0.0722*B; % 标准亮度计算
% H =  -0.11457*R - 0.38543*G - 0.5*B;      % 色度H通道
% M = 0.5*R - 0.45414*G - 0.04586*B;      % 色度M通道

end


%% 增强颜色相似性计算
function colorSim = enhanced_color_similarity(H1, M1, H2, M2, C3)

MSIM=abs((2 * (M1) .* (M2) +C3)./((M1).^2 + (M2).^2 +C3));
NSIM=abs((2 * (H1) .* (H2) +C3)./((H1).^2 + (H2).^2 + C3));
directionSim=MSIM.*NSIM;
colorSim=directionSim;
end

%% 改进的下采样函数
function [L_ds, H_ds, M_ds] = mdsi_downsample(L, H, M)
% 输入为预分离的L/H/M通道
[rows, cols] = size(L);
minDim = min(rows, cols);
F = max(1, round(minDim/256));

avgKernel = fspecial('average', F);
L_ds = imfilter(L, avgKernel, 'same');
H_ds = imfilter(H, avgKernel, 'same');
M_ds = imfilter(M, avgKernel, 'same');

L_ds = L_ds(1:F:end, 1:F:end);
H_ds = H_ds(1:F:end, 1:F:end);
M_ds = M_ds(1:F:end, 1:F:end);
end

%% 画面细节累积（仅亮度通道）
% 输入：亮度矩阵L（范围[0,100]）
% 输出：优化后的显著性图[0,1]
function saliencyMap = visualSaliency(L)
% 固定参数配置（直接内置）
 resizeFactors = [1.0, 0.5];   % 双尺度处理
% resizeFactors = 1.0;   % 单尺度处理
gradientWeights = 0.7; % 梯度权重
gaussianSigmas = [1.5, 3];    % 高斯滤波参数

%% 多尺度特征提取
saliencyMaps = cell(1,2);
for scale = 1:2
    currentL = imresize(L, resizeFactors(scale));
    [h, w] = size(currentL);

    %% 频谱残差特征（优化窗口大小）
    fftImg = fft2(single(currentL));
    logAmp = log(abs(fftImg) + eps);
    
    % 自适应高斯窗口
    winSize = 5 + 2*(scale-1); % 尺度越大窗口越大
    gKernel = fspecial('gaussian', winSize, winSize/3);
    avgLogAmp = conv2(logAmp, gKernel, 'same');
    
    spectralFeature = abs(ifft2(exp((logAmp - avgLogAmp) + 1i*angle(fftImg)))).^2;
    
    %% 梯度特征（Sobel算子+方向优化）
    [Gx, Gy] = imgradientxy(currentL, 'sobel');
    Gmag = sqrt(Gx.^2 + Gy.^2);
    
    % 四方向权重分配
    theta = atan2d(Gy, Gx);
    dirWeight = 0.5*ones(size(theta));
    dirWeight(theta >= -45 & theta < 45) = 1.0;   % 水平
    dirWeight(theta >= 135 | theta < -135) = 0.8; % 垂直
    gradientFeature = Gmag .* dirWeight;%增强水平和垂直边缘响应。
    
    %% 局部对比度特征（动态窗口）
    contrastFeature = stdfilt(currentL, true(7 + 2*(scale-1)));
    
    %% 特征融合
    maxGrad = max(gradientFeature(:)) + eps;
    maxContrast = max(contrastFeature(:)) + eps;
    saliency = zeros(h, w, 'single');
    saliency(:) = spectralFeature*gradientWeights + ...
              (gradientFeature/maxGrad)*0.6 + ...
              (contrastFeature/maxContrast)*0.4;
    saliencyMaps{scale} = imresize(saliency, size(L));
end

%% 多尺度融合
saliencyMap = 0.7*saliencyMaps{1} + 0.3*saliencyMaps{2};

%% 后处理优化
% 自适应高斯混合
sigmaBase = min(size(L))/200;
saliencyMap = 0.6*imgaussfilt(saliencyMap, gaussianSigmas(1)+sigmaBase) + ...
              0.4*imgaussfilt(saliencyMap, gaussianSigmas(2)+sigmaBase);%通过两次高斯滤波混合增强平滑性
%% 动态Sigmoid增强
mu = mean(saliencyMap(:));
sigma = std(saliencyMap(:));
saliencyMap = 1./(1 + exp(-8/(sigma+eps)*(saliencyMap - mu)));%使得高对比度区域（σ小）增强更陡峭，低对比度（σ大）增强平缓。

 end