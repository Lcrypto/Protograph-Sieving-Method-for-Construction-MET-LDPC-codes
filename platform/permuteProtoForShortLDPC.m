function [proto] = permuteProtoForShortLDPC(proto, defPun, g, iterationInPexit, targetBerPexit, layered,Samples)
[m, n] = size(proto);
infoCols = n -   m;

n1 = infoCols * 3 + defPun;
m1 = (n1 - infoCols) ;

smallProto = proto(end - m1 + 1:end, end - n1 + 1:end);
freeNodes = n1 - infoCols + 1:n1 - defPun;
punNodes = [n1 - defPun + 1:n1];
shortNodes = [];
punIt = [];%1:3:100, 3:3:100];
for it = 1:length(freeNodes)-1
    bestThreshold = 100;
    bestNodeId = -1;
    
    for nodeId = freeNodes
        if sum(it == punIt) > 0
            threshold = getThresholdsRobust(iterationInPexit, smallProto, -10, 10, ...
                {[punNodes, nodeId]}, {shortNodes},1,1, targetBerPexit, g, layered,Samples);
        else
            threshold = getThresholdsRobust(iterationInPexit, smallProto, -10, 10, ...
                {punNodes}, {[shortNodes, nodeId]},1,1, targetBerPexit, g, layered,Samples);
        end
        if threshold < bestThreshold
            bestThreshold = threshold;
            bestNodeId = nodeId;
        end
    end
    if sum(it == punIt) > 0
        punNodes = [punNodes, bestNodeId];
    else
        shortNodes = [shortNodes, bestNodeId];
    end        
    freeNodes = freeNodes(freeNodes ~= bestNodeId);
end

perm = [shortNodes, freeNodes, punNodes(end:-1:1)]
proto(:, end - infoCols + 1:end) = proto(:, perm + n - n1);
end