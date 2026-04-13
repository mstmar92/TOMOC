function [Child] = Crossover(Parent1, Parent2,...
                             N, Child)

Pc = 0.8;
for i = 1 : N
    if(rand <= Pc)
        Child.Chromosome(i) = Parent1.CommsNodes.Chromosome(i);

    else
        Child.Chromosome(i) = Parent2.CommsNodes.Chromosome(i);

    end;
end;
