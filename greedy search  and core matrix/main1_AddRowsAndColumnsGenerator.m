clear all
%%___no need to change_____________________________________________________
delete(gcp('nocreate'))
layered = 0;
%%=========================================================================

myStream=RandStream('mlfg6331_64');
rand(myStream);
degreeOfParallel = 6;
protoVersion = 0;
approximation = 1;
iterationInPexit = 50;
targetBerPexit = 1e-8;
penForDisbCoef = 0.01; 

%% proto input and settings

%%___use this to load a protograph from file_______________________________
% protofile = 'ProtoMatrices\testproto.txt';
% fproto = fopen(protofile, 'r');
% proto = [];
% while feof(fproto)==0
%     line = strsplit(fgetl(fproto), '\t');
%     protoline = [];
%     for i = 1:length(line)
%         protoline = [protoline, str2double(line{i})];
%     end
%     proto = [proto; protoline];
% end
% rcore = 0;
% ccore = 0;
%%=========================================================================

%%___use this to add 'naddcol' columns to the 'core'_______________________
core = [
0	0	0	0	0	0	1	1
1	0	0	0	1	0	1	0
1	1	1	1	0	1	0	0
1	0	1	1	0	0	0	1
1	0	1	0	1	1	0	0
0	0	0	0	0	0	1	0
];
nadd = 2;
[rcore, ccore] = size(core);
shitColsInCore = []; %[2 4 5 6 8];
proto = [core, ones(rcore,nadd)]; % use it to add columns to the core
% proto = [core, ones(nadd, ccore)]; % use it to add rows to the core
%%=========================================================================

proto = proto > 0;
[rows, cols] = size(proto);
rowStartBorder = 1; %rcore + 
rowEndBorder = rows; %
colStartBorder = ccore+1; %cols - (cols - rows) + 1; 
colEndBorder = cols; %cols;
%codeRates = [1/3, 1/4, 1/5, 1/6];
%rowsCollection = [22 32 42 52 62];
codeRates = [1/2]; % target code rate
rowsCollection = rows; %[11];
rowsForPexit = rows; %[11];
punByRates = [2];
shortByRates = [0];
coefFor4cyclesByRates = [0.01] * 0;
coefForOnesByRates = [0.001] * 666;
numPunctNodes = zeros(1, length(rowsForPexit));
numShortNodes = zeros(1, length(rowsForPexit));
coefForOnes = zeros(1, length(rowsForPexit));
coefFor4cycles = zeros(1, length(rowsForPexit));
columnsForPexit = cols + (rowsForPexit - rows) * (1 + 0);
columnsCollection = cols + (rowsCollection - rows) * (1 + 0);
for j = 1:length(rowsForPexit)
    maskMatrix{j} = zeros(rowsForPexit(j), columnsForPexit(j));
    maskMatrix{j}(rowStartBorder:rowEndBorder,colStartBorder:colEndBorder) = 1;
%     maskMatrix{j}(1:rowEndBorder,shitColsInCore) = 1;
end
%parpool(degreeOfParallel);
for i=1:length(rowsForPexit)
    for j=length(rowsCollection):-1:1
        if (rowsForPexit(i) <= rowsCollection(j))
            numPunctNodes(i) = punByRates(j);
            numShortNodes(i) = shortByRates(j);
            coefForOnes(i) = coefForOnesByRates(j);
            coefFor4cycles(i) = coefFor4cyclesByRates(j);
        end
    end
end

%% check Rates
for i = 1:length(codeRates)
    assert(codeRates(i) == (columnsCollection(i) - (1+0)*rowsCollection(i) - shortByRates(i)) /...
        (columnsCollection(i) - punByRates(i) - shortByRates(i)), 'CodeRate mismatch');
end
%%
coeffCollection = ones(1, length(rowsForPexit));

%% Parameters In Procedure
maxDifProtoInParforPexit = 2 * degreeOfParallel;
failsToStopInPexit = 40;
numberOfAttemptsForOneStep = 5;

meanNumOnes = 1.0;
for j = 1:length(maskMatrix)
    noisePr(j) = meanNumOnes / (sum(sum(maskMatrix{j})));
end


startProto = proto;
indexRateDone = 0;

bestCurLiftedMatrix = {};
while indexRateDone < length(rowsForPexit)
    coeff = coeffCollection(indexRateDone + 1);
    bestObjFun = 1000;
    for attId = 1:numberOfAttemptsForOneStep     
        threshold_search_gap_max = 10;        
        if indexRateDone~= 0
            startProto = zeros(size(maskMatrix{indexRateDone + 1}));
            for i = 1:(rowsForPexit(indexRateDone + 1) - rowsForPexit(indexRateDone))
                if 0 == 1
                    startProto(i,2*(i-1) +1:2*i) = [1 1];
                else
                    startProto(i,i) = 1;
                end
            end
            startProto(end - rowsForPexit(indexRateDone) + 1:end, end - columnsForPexit(indexRateDone)+1:end) = bestProtoOnPreviousStage;
        end        
        [rows,cols] = size(startProto);
        noiseMtr = zeros(size(maskMatrix{indexRateDone + 1}));
        while (sum(sum(noiseMtr)) == 0)
            noiseMtr = (rand(rows, cols) < (noisePr(indexRateDone + 1) *2).* maskMatrix{indexRateDone + 1});
            noiseMtr = noiseMtr.* max(startProto, ones(size(startProto)));
        end
        proto = mod(startProto + noiseMtr, 2);
        infoNodes = setdiff(1:cols, 1:(1 + 0)*rows);
        punctNodesInWhile = {};
        shortNodesInWhile = {};
        for i = 1:1
            punctNodesInWhile{i} = infoNodes(end -numPunctNodes(indexRateDone +i) +1:end);
            localSet = setdiff(infoNodes, punctNodesInWhile{i});
            shortNodesInWhile{i} = localSet(1:numShortNodes(indexRateDone + i));
        end        
        
        [~, initialThreshold, ~] = thre_search(proto, [cols-1, cols], iterationInPexit);
        
        penForDisbalance = 0;
        
        globalThreshold = initialThreshold;
        globalObjFun = sum(initialThreshold .* coeff, 2) + coefForOnes(indexRateDone + 1) * sum(sum(proto))...
            + coefFor4cycles(indexRateDone + 1) * get4cyclesByProto(proto > 0) + penForDisbalance * penForDisbCoef;
        threshold_search_gap = threshold_search_gap_max;
        noImpr = 0;
        tempr = 1;
        allIterations = 1;
        localBestProto = proto
        while true
            allIterations = allIterations + 1;
            tempr = 0.1 / allIterations;
            new_proto = zeros(rows, cols, 0);
            while size(new_proto, 3) < maxDifProtoInParforPexit
                noiseMtr = zeros(rows, cols);
                while (sum(sum(noiseMtr)) == 0)
                    noiseMtr = (rand(rows, cols) < (noisePr(indexRateDone + 1) * (1 + noImpr / (0.4 * failsToStopInPexit+1))).* maskMatrix{indexRateDone + 1});
                                        
                end
                newMtr = mod(proto + noiseMtr, 2);
                if sum(sum(newMtr, 1)==0)
                    continue;
                end
                new_proto(:, :, 1 + size(new_proto, 3)) = newMtr;
            end
            thresholds = 100* ones(1, size(new_proto,3));
             for i=1:size(new_proto, 3)
           % parfor i=1:size(new_proto, 3)
                [~, thresholds(:, i), ~] = thre_search(new_proto(:,:,i), [cols-1, cols], iterationInPexit)
                
                objFun(i) = sum(thresholds(:, i).*coeff, 2) + coefForOnes(indexRateDone + 1) * sum(sum(new_proto(:,:,i))) + ...
                    coefFor4cycles(indexRateDone + 1) * get4cyclesByProto(new_proto(:,:,i) > 0) + penForDisbalance * penForDisbCoef;
            end
            [localObjFun,ind] = min(objFun);
            if (abs(localObjFun-globalObjFun)<1e-4 || localObjFun>globalObjFun) && (noImpr == failsToStopInPexit)
                break;
            end
            threshold_search_gap = min(threshold_search_gap_max, max(abs(localObjFun-globalObjFun)*2, 0.01))
            if (globalObjFun > 20)
                threshold_search_gap = 5;
            end
            if (globalObjFun > localObjFun)
                globalObjFun = localObjFun;
                globalThreshold = thresholds(:, ind);
                proto = new_proto(:,:,ind);
                localBestProto = proto;
                
                noImpr = 0;
            else
                p = rand();
                if (localObjFun < 30) && (p < exp(100 / (1 + 0.1*noImpr) *(globalObjFun - 1e-4 - localObjFun) / tempr))
                    
                    proto = new_proto(:,:,ind);
                    noImpr = 0;
                else
                    if mod(noImpr, 10) == 0
                        
                    end
                    noImpr = noImpr + 1;
                end
            end
        end
        
        
        fprintf('rows = %d cols = %d threshold = %f\n', rowsForPexit(indexRateDone + 1), columnsForPexit(indexRateDone + 1), globalThreshold);
        
        if globalObjFun < bestObjFun
            bestObjFun = globalObjFun;
            thresholdInBestObjFun = globalThreshold
            bestProto = proto;
        end
        
    end
    bestProtoOnPreviousStage = bestProto;
    fprintf('rows = %d cols = %d best object function = %f threshold = %f\n', ...
        rowsForPexit(indexRateDone + 1), columnsForPexit(indexRateDone + 1), bestObjFun, thresholdInBestObjFun);
    indexRateDone = indexRateDone + 1;
end
%%
% [bestProto] = permuteProtoForShortLDPC(bestProto, punByRates(1), 0, iterationInPexit, targetBerPexit, layered);

[rows, cols] = size(bestProto);
fid = fopen(['proto', num2str(cols), 'x', num2str(rows), 'ver', num2str(version), '.txt'], 'w');
fprintf(fid, '%d\t%d\n', cols, rows);
for r = 1:rows
    fprintf(fid, '%d\t', bestProto(r, :));
    fprintf(fid, '\n');
end
fclose(fid);


