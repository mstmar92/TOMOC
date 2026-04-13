function [Results]= MOEAD_Overlapping(A, N, m, ...
                          IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
                          NumberOfProteinsInComplexes, NumberOfKnownProteinsInComplexes, ...
                          Population, PopulationSize,...
                          Heuristic, Pm, ...
                          KnownProteins, ProteinLabelPairs)

                                   
                                   
%-------------------------------------------------------------------------%
%                        MOEA/D Parameter Settings                        %
%-------------------------------------------------------------------------%
ObjectiveDimension = 2; % Number of objectives to be optimized

ParameterDimension = length(A); % similar vars.

% All models seeks for minimization                                               
IdealPoint(1:ObjectiveDimension) = 10000000; % all other MOO models, including our model, seek for Minimization


IndivPoint = [];
 
Params.PopSize = PopulationSize;
 
Params.Neighbour = 5;
 
Params.Generations = 100; %500;%250 
 
Params.DMethod = 'ts'; % there are other varients of decomposition
 
SubProblem = InitWeights(Params.PopSize, Params.Neighbour, ObjectiveDimension);
 
Params.PopSize = length(SubProblem);

Params.Mut.MinSize = 3;

Params.Mut.pAddTriangle = 0.6;

Params.Mut.pPrune = 0.4;

Params.Mut.MaxComplexAddsPerNode = 2;

Params.Mut.MaxMovesPerNode = 1;


%-------------------------------------------------------------------------%
%                                  MOEA/D                                 %
%-------------------------------------------------------------------------%

% Update Reference points
[IdealPoint, IndivPoint] = UpdateReferenceOverlapping(IdealPoint, IndivPoint, ...
                                           Population, Params.PopSize);
                                            

% Initialize the external population, which finally holds the non-dominated solutions
Results = [];
NearParetoOptimalSet = [];

[NearParetoOptimalSet] = UpdateNearParetoOptimalSetOverlapping(NearParetoOptimalSet, ...
                                                    Population, Params.PopSize, ...
                                                    ObjectiveDimension);
                                                    
                                                 
Results(1).NearParetoOptimalSet = NearParetoOptimalSet;                                                 

% MOEA/D Loop
for GenerationCounter = 2 : Params.Generations
    GenerationCounter
    IdealPoint
    Child = Population;
    for ProblemCounter = 1 : Params.PopSize
        % 5.1 Selection Operator
        rand1 = floor(rand*(Params.Neighbour-1)+1);
        rand2 = floor(rand*(Params.Neighbour-1)+1);
        
        P1 = SubProblem(ProblemCounter).Neighbour (rand1);
        P2 = SubProblem(ProblemCounter).Neighbour (rand2);
        
        Parent1(ProblemCounter) = Population(P1);
        Parent2(ProblemCounter) = Population(P2);
    end;
        % 5.2 Crossover Operator
        for ProblemCounter1 = 1 : Params.PopSize
        [Child(ProblemCounter1)] = Crossover(Parent1(ProblemCounter1), Parent2(ProblemCounter1), ...
                                            N, ...
                                            Child(ProblemCounter1));
        end;

        % 5.3 SOR Mutation Operator
    for ProblemCounter3 = 1 : Params.PopSize
                                        
        Child(ProblemCounter3) = Mutation_SOR_Overlapping(Child(ProblemCounter3), A, ...
                                         IndicesInteractionProtein, NumInteractionProtein, ...
                                         Pm, Params.Mut);  

    end;
                    

    [Child] = Individual2CmplxDecodingOverlapping(ProteinLabelPairs, ...
                                            Child, PopulationSize);
    [Child] = ComputeFitnessOverlapping(A, N, m, ...
                                       IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
                                       Child, PopulationSize);
    
    [IdealPoint, IndivPoint] = UpdateReferenceOverlapping(IdealPoint, IndivPoint, ...
                                                Child, Params.PopSize);
                                                
                                            
    [Population] = UpdateProblemOverlapping(Population, Child, ...
                                  SubProblem, ...
                                  Params, ...
                                  IndivPoint, ...
                                  ObjectiveDimension);
                                 
    [NearParetoOptimalSet] = UpdateNearParetoOptimalSetOverlapping(NearParetoOptimalSet, ...
                                                         Child, Params.PopSize, ...
                                                         ObjectiveDimension);
                                                         
                                           
                                                 
    Results(GenerationCounter).NearParetoOptimalSet = NearParetoOptimalSet;                                                 
   
end;                                   