function [NearParetoOptimalSet] = Dominate(Child, ...
                                           NearParetoOptimalSet, nd, ...
                                           ObjectiveDimension)

if ObjectiveDimension == 2
    if ( (Child.Obj(1) <= NearParetoOptimalSet(nd).Obj(1)) && ...
         (Child.Obj(2) <= NearParetoOptimalSet(nd).Obj(2)) && ...
         ( (Child.Obj(1) < NearParetoOptimalSet(nd).Obj(1)) || ...
           (Child.Obj(2) < NearParetoOptimalSet(nd).Obj(2)) ) )

        NearParetoOptimalSet(nd).Dominated = 'T';
    else
        NearParetoOptimalSet(nd).Dominated = 'F';
    end
end

end
