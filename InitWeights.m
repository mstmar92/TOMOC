function SubProblem = InitWeights (PopSize, Neighbour, ObjectiveDimension)
%--------------------------------------------------------------------------
% init_weights function initialize a population of subproblems structure
% with the generated decomposition weight and the neighbourhood relationship.

SubProblem = [];
for WeightVectorCounteri = 1 : PopSize    
    P = struct ('Weight', [], 'Neighbour', [], 'Optimal', Inf, 'OptPoint', [], 'Curpoint', []);
    if (ObjectiveDimension == 2)
        weight(1) = WeightVectorCounteri / PopSize;
        weight(2) = (PopSize-WeightVectorCounteri) / PopSize;
        P.Weight = weight;
        SubProblem = [SubProblem P];
    else
        if (ObjectiveDimension == 3)
            for WeightVectorCounterj = 1:PopSize
                if(WeightVectorCounteri + WeightVectorCounterj <= PopSize)
                   WeightVectorCounterk = PopSize - WeightVectorCounteri - WeightVectorCounterj; 
                   weight(1) = WeightVectorCounteri / PopSize;
                   weight(2) = WeightVectorCounterj / PopSize;
                   weight(3) = WeightVectorCounterk / PopSize;
                   P.Weight = weight;
                   SubProblem = [SubProblem P];
                end;
            end;
        end;
    end;   
end;


Length = length (SubProblem);
%Length = PopSize;
DistanceMatrix = zeros(Length, Length);
for i = 1: Length
    
    for j= i+1:Length
        
        A = SubProblem(i).Weight;
        B= SubProblem(j).Weight;
        DistanceMatrix(i,j) = (A-B)*(A-B)';
        DistanceMatrix(j,i) =  DistanceMatrix(i,j);
        
    end;
    
    [s,sindex] = sort ( DistanceMatrix(i,:));
    
    SubProblem(i).Neighbour = sindex(1:Neighbour)';
    
end;
SubProblem   