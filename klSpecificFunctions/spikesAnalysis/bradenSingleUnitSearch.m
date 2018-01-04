function bradenSingleUnitSearch(spks,varargin)

close all;

global Target_ SRT Correct_ SAT_

flankTarg = 0;
flankDist = 1;
doAll = 1;
satCrit = [1,2];

% Decode varargin here?
varStrInd = find(cellfun('isclass',varargin,'char'));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'ft'}
            flankTarg = varargin{varStrInd(iv)+1};
        case {'fd'}
            flankDist = varargin{varStrInd(iv)+1};
        case {'a','-a'}
            doAll = varargin{varStrInd(iv)+1};
    end
end

if flankTarg
    bsxTarg = -1:1;
else
    bsxTarg = 0;
end
if flankDist
    bsxDist = -1:1;
else
    bsxDist = 0;
end

% Get SDFs
vCheck = 50:150; mCheck = -50:0;
[vSDF,vTimes] = klSpkRatev2(spks-repmat(Target_(:,1),1,size(spks,2)),'-q',1);
[mSDF,mTimes] = klSpkRatev2(spks-repmat(Target_(:,1)+SRT(:,1),1,size(spks,2)));
myTimes = find(vTimes >= 50 & vTimes <= nanmean(SRT(:,1))-50);

% Convert target vector to x and y
spX = [1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2,0,sqrt(2)/2];
spY = [0,sqrt(2)/2,1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2];

xLoc = nan(size(vSDF,1),1);
yLoc = nan(size(vSDF,1),1);
for il = 1:length(spX)
    xLoc(Target_(:,2)==(il-1)) = spX(il);
    yLoc(Target_(:,2)==(il-1)) = spY(il);
end
genCrit = Correct_(:,2)==1;
if ~isempty(SAT_)
    genCrit = genCrit & ismember(SAT_(:,1),satCrit);
end


% Get average responses per location
uLocs = nunique(Target_(:,2));
vResp = nan(1,length(uLocs));
mResp = nan(1,length(uLocs));
for il = 1:length(uLocs)
    vResp(il) = nanmean(nanmean(vSDF(Correct_(:,2)==1 & Target_(:,2)==uLocs(il),myTimes),1));
    mResp(il) = nanmean(nanmean(mSDF(Correct_(:,2)==1 & Target_(:,2)==uLocs(il),myTimes),1));
end
if ~any(~isnan(vResp)), vrf = nan; else vrf=uLocs(vResp==max(vResp)); end; if isnan(vrf) || isempty(vrf), badV = 1; end
if ~any(~isnan(mResp)), mrf = nan; else mrf=uLocs(mResp==max(mResp)); end; if isnan(mrf) || isempty(mrf), badM = 1; end

figure();
spInds = [6,3,2,1,4,7,8,9];
for i = 1:8
    sp(i)=subplot(3,3,spInds(i));
    plot(vTimes,nanmean(vSDF(genCrit & Target_(:,2)==(i-1),:),1));
    set(gca,'XLim',[-50,nanmean(SRT(:,1))-50],'box','off');
end
linkaxes(sp,'y');
midAx = subplot(3,3,5);
axLoc = get(midAx,'Position');
delete(midAx);
polaraxes('Position',axLoc);
polar(klDeg2Rad(45.*[uLocs;uLocs(1)]),[vResp,vResp(1)]');
% rticks([]); thetaticks([]);
xticks([]); yticks([]);

% Get respX and respY
if doAll
    allVResp = nanmean(vSDF(:,myTimes),2);
    respX = xLoc.*allVResp; respY = yLoc.*allVResp;
    mnX = nanmean(respX(genCrit));
    mnY = nanmean(respY(genCrit));
    for il = 1:length(uLocs)
        vResp(il) = nansum(allVResp(genCrit & Target_(:,2)==uLocs(il)))./sum(genCrit);
    end
else
    respX = spX.*[vResp];
    respY = spY.*[vResp];
    respR = sqrt(respX.^2+respY.^2);
    respTh = klDeg2Rad(45.*[uLocs]);
    mnX = nanmean(respX);
    mnY = nanmean(respY);
end

mnR = sqrt(mnX^2+mnY^2);
mnTh = klRad2Deg(abs(atan(mnY/mnX)));
if (mnX >= 0) && (mnY>=0)
elseif (mnX < 0) && (mnY >= 0)
    mnTh = 180-mnTh;
elseif (mnX < 0) && (mnY < 0)
    minTh = 180+mnTh;
else
    mnTh = 360-mnTh;
end

figure();
polAx = polaraxes;
polar(klDeg2Rad([uLocs;uLocs(1)].*45),[vResp,vResp(1)]'./max(vResp)); hold on;
polar([0,klDeg2Rad(mnTh)],[0,mnR]);


% Get running average SDF
avSDF = klRunningAvv2(vSDF,5);
myTimes = find(vTimes >= -50 & vTimes <= nanmean(SRT(:,1))-50);

targIn = avSDF(genCrit & ismember(Target_(:,2),bsxfun(@plus,vrf,bsxTarg)),myTimes);
targOut = avSDF(genCrit & ~ismember(Target_(:,2),bsxfun(@plus,vrf,bsxDist)),myTimes);

% % Calculate AUC
subTimes = vTimes(myTimes(1:5:end));
auc = nan(1,length(subTimes));
cis = nan(length(subTimes),2);
for ic = 1:length(subTimes)
    [auc(ic),~,~,cis(ic,:)] = klROCv1(targIn(:,1+(5*(ic-1))),targOut(:,1+(5*(ic-1))));
end
% 
% figure();
% subplot(3,1,[1,2]); hold on;
% plot(vTimes(myTimes),nanmean(targIn,1),'k','linewidth',3);
% plot(vTimes(myTimes),nanmean(targOut,1),'k','linewidth',1);
% subplot(3,1,3); hold on;
% plot(subTimes,auc);
% plot(subTimes,cis(:,1),'k','linestyle','--');
% plot(subTimes,cis(:,2),'k','linestyle','--');
    
    keyboard
    