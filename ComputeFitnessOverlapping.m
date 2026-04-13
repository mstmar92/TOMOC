function [Population] = ComputeFitnessOverlapping(A, N, m, ...
                                                 IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
                                                 Population, PopulationSize)
                                             
         

 
for i = 1:PopulationSize
Population(i).Obj = [];
Population(i).Obj1_Conductance = 0;
Population(i).Obj2_TriLack = 0;
Population(i).Kv = 0;
Population(i).avgTri = 0;
end

for IndividualCounter = 1:PopulationSize
Population(IndividualCounter) = ComputeOurFitnessCollectionOverlapping( ...
        A, N, IndicesInteractionProtein, NumInteractionProtein, ...
        Population(IndividualCounter));
end