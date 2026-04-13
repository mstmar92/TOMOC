function [NearParetoOptimalSet] = EraseFromNearParetoOptimalSetOverlapping(NearParetoOptimalSet)

NearParetoOptimalSetTemp = [];
cnt = 0;

for i = 1:length(NearParetoOptimalSet)
    if (NearParetoOptimalSet(i).Dominated ~= 'T')
        cnt = cnt + 1;
        NearParetoOptimalSetTemp(cnt).Chromosome = NearParetoOptimalSet(i).Chromosome;
        NearParetoOptimalSetTemp(cnt).Dominated  = NearParetoOptimalSet(i).Dominated;
        NearParetoOptimalSetTemp(cnt).CmplxID    = NearParetoOptimalSet(i).CmplxID;

        NearParetoOptimalSetTemp(cnt).Obj        = NearParetoOptimalSet(i).Obj;
        NearParetoOptimalSetTemp(cnt).Obj1_Conductance = NearParetoOptimalSet(i).Obj1_Conductance;
        NearParetoOptimalSetTemp(cnt).Obj2_TriLack     = NearParetoOptimalSet(i).Obj2_TriLack;
    end
end

NearParetoOptimalSet = NearParetoOptimalSetTemp;
end
