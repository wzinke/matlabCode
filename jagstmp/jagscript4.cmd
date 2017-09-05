load dic
model in "/home/loweka/git/matlabCode/tdtExtraction/Rate_1.txt"
data in jagsdata.R
compile, nchains(1)
parameters in jagsinit4.R
initialize
update 1000
monitor set theta, thin(1)
monitor deviance
update 5000
coda *, stem('CODA4')
