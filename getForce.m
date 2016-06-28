% written by: Stu Blair
% date: Jan 16, 2012
% purpose: get force vector for all solid nodes in an LBM domain

function F = getForce(solidNodes,streamTgtMat,LatticeSpeeds,bb_spd,fIn)

% determine number of solid nodes
numSolidNodes=length(solidNodes);

%deterimne number of dimensions
[numDim,numSpd]=size(LatticeSpeeds);

F = zeros(numSolidNodes,numDim);
% get forceContribMat - map of force contributors for each solid node
forceContribMat = zeros(numSolidNodes,numSpd);
for spd = 2:numSpd % stationary speed never contributes
    forceContribMat(:,spd)=~ismember(streamTgtMat(solidNodes,spd),solidNodes);
end

e_alpha = LatticeSpeeds';

for spd=2:numSpd
    f1 = fIn(solidNodes,spd)+fIn(streamTgtMat(solidNodes,spd),bb_spd(spd));
    f2 = (f1).*forceContribMat(:,spd);
    f3 = kron(f2,e_alpha(bb_spd(spd),:));
    F = F + f3;
  
end
