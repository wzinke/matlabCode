%% This function implements Landmark-based spectral clustering (LSC) as defined by Chen & Cai (2011) (via Nguyen et al 2014)

function [out] = klLSCv1(inVals,varargin)

% Set defaults
p = min([1000,size(inVals,1)/2]); % Set number of landmarks
pType = 'rand'; % How to get landmarks
r = 5; % # nearest neighbors
gSD = 5; % kernel bandwidth
k = 1:6;

%% Step 1: create landmarks p
switch pType
    case 'kmeans'
        [~,landmarks] = kmeans(inVals,p);
    case 'rand'
        randSet = randperm(size(inVals,1));
        landmarks = inVals(randSet(1:p),:);
end

% Get a distance matrix just in case
landDists = EuDist2(inVals,landmarks);

%% Step 2: produce sparce affinity matrix Z (dim=pxn) between data points
% and landmark points with the affinity computed as:
% 
% z(j,i) = [Kh(xi,uj)]/[sum(j')(Kh(xi,uj'))
% 
% where Kh() is a kernel function with a bandwiths h and j, where j is in
% Ui which is, in turn, a submatrix of U composed of the r nearest
% landmarks of xi. uj are landmark points
%
% Usually a Gaussian kernel is used:
%    Kh(xi,uj) = exp(-||xi-uj||^2/2h^2)

% Set up kernel
gKern = @(xi,uj,h) exp(-1*(nansum((xi-uj).^2,2))/(2*h^2));

% Loop through to fill Z (hopefully I can think of a way to remove the loop
% later?
Z = zeros(p,size(inVals,1));
for i = 1:size(inVals,1),
    % Get r nearest neighbors
    [~,sortInds] = sort(landDists(i,:),'ascend');
    rNearestLandmarks = landmarks(sortInds(1:r),:);
    for j = 1:p,
        Z(j,i) = gKern(inVals(i,:),landmarks(j,:),gSD)/sum(gKern(repmat(inVals(i,:),r,1),rNearestLandmarks,gSD));
    end
end

%%Step 3: Calculate the first k eigenvectors of ZZt, denoted by A =
%%[a1,a2,...,ak]

ZZT = Z*Z';
[vals,diagMat] = eig(ZZT);

A = diagMat(:,1:max(k));

%% Step 4: Calculate B=[b1,b2,...,bk] (an nxp matrix) by B' = sum(-1A'*Zhat) where Zhat = D^(-1/2)*Z and D is row sum of Z
D = nansum(Z,1);
Zhat = (D.^(-1/2))*Z;


