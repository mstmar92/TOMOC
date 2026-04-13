function [NearParetoOptimalSet] = UpdateNearParetoOptimalSetOverlapping(NearParetoOptimalSet, ...
                                                              Population, PopSize, ...
                                                              ObjectiveDimension)

for i = 1:PopSize
    Flag = 1;

    for nd = 1:length(NearParetoOptimalSet)

        if isequal(NearParetoOptimalSet(nd).CmplxID, Population(i).CommsNodes.CommsNodes)
            NearParetoOptimalSet(nd).Dominated = '';
            Flag = 0;
            break;
        end

        if ObjectiveDimension == 2
            if ( (NearParetoOptimalSet(nd).Obj(1) <= Population(i).Obj(1)) && ...
                 (NearParetoOptimalSet(nd).Obj(2) <= Population(i).Obj(2)) && ...
                 ( (NearParetoOptimalSet(nd).Obj(1) < Population(i).Obj(1)) || ...
                   (NearParetoOptimalSet(nd).Obj(2) < Population(i).Obj(2)) ) )
                NearParetoOptimalSet(nd).Dominated = '';
                Flag = 0;
                break;
            end
        end

        NearParetoOptimalSet = Dominate(Population(i), NearParetoOptimalSet, nd, ObjectiveDimension);
    end

    if Flag
        NearParetoOptimalSet = EraseFromNearParetoOptimalSetOverlapping(NearParetoOptimalSet);

        L = length(NearParetoOptimalSet);
        NearParetoOptimalSet(L+1).Chromosome = Population(i).Chromosome;
        NearParetoOptimalSet(L+1).Dominated  = '';
        NearParetoOptimalSet(L+1).CmplxID    = Population(i).CommsNodes.CommsNodes;

        NearParetoOptimalSet(L+1).Obj        = Population(i).Obj;   % <-- new
        NearParetoOptimalSet(L+1).Obj1_Conductance = Population(i).Obj(1);
        NearParetoOptimalSet(L+1).Obj2_TriLack     = Population(i).Obj(2);
    end
end

end
