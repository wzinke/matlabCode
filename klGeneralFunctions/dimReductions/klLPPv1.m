%% This function consolidates the LPP and constructW and option creation
%  from Cai et al/He et al's LPP code which is called below. This allows
%  the user to start directly with the input data without fussing with
%  these things necessarily

function [outVals, eigVector, eigValue] = klLPPv1(inData,varargin)


% Set defaults
options.Metric = 'Euclidean';
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'HeatKernel';
options.t = 5;
options.PCARatio = .99;

% Decode varargin


% Set up W (affinity matrix)
W=constructW(inData,options);

% Now call LPP
[eigVector,eigValue] = LPP(W, options, inData);

% Do the transformation
outVals = inData*eigVector;