function [diffPerc,clustOverlap] = klContingDev(respClusts,waveClusts),


u1 = unique(respClusts);
if iscell(u1(1)),
    if ischar(u1{1}),
        u1(strcmp(u1,'')) = [];
        tmp1 = nan(size(respClusts));
        for i1 = 1:length(u1),
            tmp1(strcmp(respClusts,u1{i1})) = i1;
        end
        respClusts = tmp1;
        u1 = unique(tmp1);
    elseif isnumeric(u1{1}),
        u1 = cell2mat(u1);                   
        u1(isnan(u1)) = [];
    end
end
u1(isnan(u1)) = [];

u2 = unique(waveClusts);
if iscell(u2(1)),
    if ischar(u2{1}),
        u2(strcmp(u2,'')) = [];
        tmp2 = nan(size(waveClusts));
        for i2 = 1:length(u2),
            tmp2(strcmp(waveClusts,u2{i2})) = i2;
        end
        waveClusts = tmp2;
        u2 = unique(tmp2);
    elseif isnumeric(u2{1}),
        u2 = cell2mat(u2);                   
        u2(isnan(u2)) = [];
    end
else
    u2(isnan(u2)) = [];
end

for ir = 1:length(u1),
    for ic = 1:length(u2),
        inMat(ir,ic) = sum(respClusts == u1(ir) & waveClusts == u2(ic));
    end
end

u2 = unique(waveClusts(~isnan(waveClusts)));
u1 = unique(respClusts(~isnan(respClusts)));

for iw = 1:length(u2),
    for ir = 1:length(u1),
        clustOverlap(iw,ir) = sum(respClusts==u1(ir) & waveClusts == u2(iw));
    end
end
for iw = 1:length(u2),
    for ir = 1:length(u1),
        expect(iw,ir) = (sum(respClusts==u1(ir))./sum(isfinite(respClusts)))*(sum(waveClusts == u2(iw))./sum(isfinite(waveClusts)));
    end
end
% expect = expect(:,[1,3,4]);
allOverlap = clustOverlap./sum(clustOverlap(:));
expectAll = expect.*sum(clustOverlap(:));
% clustOverlap = clustOverlap(:,[1,3,4]);
diff = clustOverlap-expectAll;
diffPerc = diff./expectAll;