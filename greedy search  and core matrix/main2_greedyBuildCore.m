clear all;
clc;

%% Settings
proto = [
0	0	0	0	0	0	1	1
1	0	0	0	1	0	1	0
1	1	1	1	0	1	0	0
1	0	1	1	0	0	0	1
1	0	1	0	1	1	0	0
0	0	0	0	0	0	1	0
];
[m, n] = size(proto);
pun = n-1:n;
iter = 50; % number of decoding iterations
Wc = zeros(1,n);
for i = 1:n
    Wc(i) = sum(proto(:,i));
end
Wr = zeros(1,m);
for i = 1:m
    Wr(i) = sum(proto(i,1:n));
end
%========================================================================

%% Main

%%___tune this block to define mask by weights___
% badWeight = 2;
% badRowWeight = 1;
% badCols = find(Wc == badWeight);
% badRows = find(Wr <= badRowWeight);
% targetWeight = 3;
%%=============================================

%%_use this to consider whole graph without mask_
badCols = 1:n;
badRows = 1:m;
%%=============================================

freepoint = 10000; % use freepoint to chose a column after which you are
                 % ready to put '1' instead of '0' at any row or keep it >n                
bestProto = proto;
bestSNR = 100; % use this if you start from unvalid graph
% [~, bestSNR, ~] = thre_search(S, pun, iter); % use this if you start from valid
                                    % graph
for i = 1:length(badCols)
    S = bestProto;     
    if i >= freepoint
        badRows = 1:m;
    end
    col = proto(:,badCols(i));
    variants = setdiff(badRows, find(col ~= 0));
    for j = 1:length(variants)
        S2 = S;
        S2(variants(j), badCols(i)) = 1;   
        [~, SNR, ~] = thre_search(S2, pun, iter);
        if SNR < bestSNR
            bestSNR = SNR;
            bestProto = S2;
        end
    end
%     save(savefile, 'bestSNR', 'bestProto');
end












