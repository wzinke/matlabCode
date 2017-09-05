function histMat = klHist2(obs,bins1,bins2)

infBins1 = [-inf,bins1(2:end),inf];
infBins2 = [-inf,bins2(2:end),inf];

histMat = nan(length(bins2),length(bins1));

for ir = 1:(length(infBins1)-1),
    for ic = 1:(length(infBins2)-1),
        in1 = obs(:,1) >= infBins1(ir) & obs(:,1) < infBins1(ir+1);
        in2 = obs(:,2) >= infBins2(ic) & obs(:,2) < infBins2(ic+1);
        histMat(ic,ir) = sum(in1 & in2);
    end
end
% keyboard
% histMat(1,1) = sum(obs(:,1) < min(bins1) & obs(:,2) < min(bins2));
% for ib2 = 1:(length(bins2)-1),
%     exp2 = obs(:,2) >= bins2(ib2) & obs(:,2) < bins2(ib2+1);
%     histMat(1,ib2+1) = sum(exp2 & obs(:,1) < min(bins1));
% end
% histMat(1,ib2+2) = sum(obs(:,1) < min(bins1) & obs(:,2) >= max(bins2));
% 
% for ib1 = 1:(length(bins1)),
%     if ib1 == length(bins1),
%         exp1 = obs(:,1) >= max(bins1);
%     else
%         exp1= obs(:,1) >= bins1(ib1) & obs(:,1) < bins1(ib1+1);
%     end
%     histMat(ib1+1,1) = sum(exp1 & obs(:,2) < min(bins2));
%     for ib2 = 1:(length(bins2)-1),
%         exp2 = obs(:,2) >= bins2(ib2) & obs(:,2) < bins2(ib2+1);
%         histMat(ib1+1,ib2+1) = sum(exp1 & exp2);
%     end
%     histMat(ib1+1,ib2+2) = sum(exp1 & obs(:,2) >= max(bins2));
% end
% keyboard
% histMat = nan(length(bins1),length(bins2));
% histMat(1,1) = sum(obs(:,1) <= bins1(1) & obs(:,2) <= bins2(1));
% for ir = 2:length(bins1)-1,
%     histMat(ir,1) = sum(obs(:,1) >= bins1(ir) & obs(:,1) < bins1(ir+1) & obs(:,2) <= bins1(1));
%     for ic = 2:length(bins2)-1,
%         histMat(ir,ic) = sum(klInBox(obs(:,1),obs(:,2),[nanmean(bins1(ir:(ir+1))),nanmean(bins2(ic:(ic+1)))],[(bins1(ir+1)-bins1(ir))/2,(bins2(ic+1)-bins2(ic))/2]));
%     end
%     histMat(ir,length(bins2)) = sum(obs(:,1) >= bins1(ir) & obs(:,1) < bins1(ir+1) & obs(:,2) > bins2(end));
% end
% histMat(end,end) = sum(obs(:,1) > bins1(end) & obs(:,2) > bins2(end));

    
    
    
    
    
    