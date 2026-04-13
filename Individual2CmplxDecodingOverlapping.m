function [Population] = Individual2CmplxDecodingOverlapping(ProteinLabelPairs, ...
                                                            Population, PopulationSize)


for IndividualCounter = 1 : PopulationSize
    Population(IndividualCounter).CommsNodes = ComputeCmplxDecodingOverlapping(Population(IndividualCounter), ProteinLabelPairs);
end;