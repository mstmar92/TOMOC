function [Individual] = UpdateProblemOverlapping(Individual, Child, ...
                                       SubProblems, ...
                                       Params, ...
                                       IndivPoint, ...
                                       ObjectiveDimension)

% Minimization + Tchebycheff

for ProblemCounter = 1:Params.PopSize
    for NeighbourCounter = 1:Params.Neighbour

        NeighbourIndex = SubProblems(ProblemCounter).Neighbour(NeighbourCounter);

        IndividualofNeighbourIndex = Individual(NeighbourIndex);
        WeightofNeighbourIndex = SubProblems(NeighbourIndex).Weight;

        F1 = ScalarFunction(IndividualofNeighbourIndex, WeightofNeighbourIndex, ...
                            IndivPoint, ObjectiveDimension);

        F2 = ScalarFunction(Child(ProblemCounter), WeightofNeighbourIndex, ...
                            IndivPoint, ObjectiveDimension);

        if (F2 < F1)
            Individual(NeighbourIndex).Chromosome = Child(ProblemCounter).Chromosome;
            Individual(NeighbourIndex).CommsNodes.CommsNodes = Child(ProblemCounter).CommsNodes.CommsNodes;
            Individual(NeighbourIndex).Obj = Child(ProblemCounter).Obj;

            % Optional: keep named fields for debugging
            Individual(NeighbourIndex).Obj1_Conductance = Child(ProblemCounter).Obj(1);
            Individual(NeighbourIndex).Obj2_TriLack     = Child(ProblemCounter).Obj(2);
        end
    end
end

end
