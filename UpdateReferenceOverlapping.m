function [IdealPoint, IndivPoint] = UpdateReferenceOverlapping(IdealPoint, IndivPoint, ...
                                                     Population, PopulationSize)
% Minimization - update ideal point z for MOEA/D

for i = 1:PopulationSize

    % Objective 1: Conductance (min)
    if Population(i).Obj(1) < IdealPoint(1)
        IdealPoint(1) = Population(i).Obj(1);
        IndivPoint(1).Chromosome = Population(i).Chromosome;
        IndivPoint(1).CmplxID = Population(i).CommsNodes.CommsNodes;
        IndivPoint(1).Ideal = IdealPoint(1);   % store ideal value
        IndivPoint(1).Obj = Population(i).Obj;
    end

    % Objective 2: Triangle Lack (min)
    if Population(i).Obj(2) < IdealPoint(2)
        IdealPoint(2) = Population(i).Obj(2);
        IndivPoint(2).Chromosome = Population(i).Chromosome;
        IndivPoint(2).CmplxID = Population(i).CommsNodes.CommsNodes;
        IndivPoint(2).Ideal = IdealPoint(2);
        IndivPoint(2).Obj = Population(i).Obj;
    end

end
end
