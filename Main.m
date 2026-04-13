% -------------------------------------------------------------------------------------------------------%
%    Multi-objective EA for Overlapping complex Detection in Protein-Protein Interaction Networks        %
% -------------------------------------------------------------------------------------------------------%

clc;
clear;
NetworkNumber = input('Enter PPINs Network: ');

if (NetworkNumber == 1)
    load('DataSets\PPI\1-Protein-Yeast-D1-Files.mat');
    load('DataSets\Complex\Complex-D1-Files.mat');
elseif (NetworkNumber == 2)
    load('DataSets\PPI\2-Protein-Yeast-D2-Files.mat');
    load('DataSets\Complex\Complex-D2-Files.mat');
elseif (NetworkNumber == 4)
    DatasetName = "CollinsCYC";
    load('DataSets\PPI\4-Collins.mat');
    load('DataSets\Complex\4-CYC2008.mat'); 
elseif (NetworkNumber == 5)
    load('DataSets\PPI\5-Collins.mat');
    load('DataSets\Complex\5-MIPS.mat');    
end

% 3) Initialization Parameters
MaxRun = 10;
Pm = 0.2;
Heuristic = 1;

m = sum(sum(A)) / 2;

PopulationSize = 100;
ResultsGroup = [];   
PopulationGroup = []; 

for RunNumber = 1 : MaxRun
    [PopulationGroup(RunNumber).Population] = CreatePopulationOverlapping(ProteinLabelPairs, ...
                                             PopulationSize);
end;

for RunNumber = 1 : MaxRun
    RunNumber
                                 
    % Step 3: Chromosome (Solution) decoding
      [PopulationGroup(RunNumber).Population] = Individual2CmplxDecodingOverlapping(ProteinLabelPairs, ...
                                            PopulationGroup(RunNumber).Population, PopulationSize);                                  
                                                 
 
    % Step 4: Compute collection of Fitness functions.
    [PopulationGroup(RunNumber).Population] = ComputeFitnessOverlapping(A, N, m, ...
                                                                       IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
                                                                       PopulationGroup(RunNumber).Population, PopulationSize);  
    % Step 5: Compute ResultsGroup by MOEAD                          
    [ResultsGroup(RunNumber).Results] = MOEAD_Overlapping(A, N, m, ...
                                              IndicesInteractionProtein, NumInteractionProtein, MaxNumInteractionProtein, ...
                                              NumberOfProteinsInComplexes, NumberOfKnownProteinsInComplexes, ...
                                              PopulationGroup(RunNumber).Population, PopulationSize, ...
                                              Heuristic, Pm, ...
                                              KnownProteins, ProteinLabelPairs);
                                 

end; % for   RunNumber 1 to MaxRuns

save(strcat(strcat('RepOverlap/MOEA/MOEAD_', 'PPI_', int2str(NetworkNumber), ...
                   '_Heuristic_', int2str(Heuristic), '_Pm_', num2str(Pm), ...
                   'NumOfRuns', int2str(MaxRun)),'.mat'), 'ResultsGroup', '-v7.3');