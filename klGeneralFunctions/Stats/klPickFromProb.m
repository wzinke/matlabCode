function outX = klPickFromProb(probs)

% First convert probs to CDF
probCum = cumsum(probs./sum(probs));

% Pick random value
randVal = rand;

% Get the last probCum <= randVal; assign randVal=0 to first stim
outX = find(probCum >= randVal,1);
if isempty(outX)
    outX = 1;
end