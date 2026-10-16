% Problem 4: step experiments, identification, then fitted-model Markov coefficients.
problemDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(problemDirectory,'..','common'),problemDirectory);
[problem04DeterministicResults, problem04DeterministicFigures] = problem04Deterministic();
saveExperiment(fullfile(problemDirectory,'results'),'deterministic_steps',problem04DeterministicResults,problem04DeterministicFigures);
[problem04StochasticResults, problem04StochasticFigures] = problem04Stochastic();
saveExperiment(fullfile(problemDirectory,'results'),'stochastic_steps',problem04StochasticResults,problem04StochasticFigures);
[problem04IdentificationResults, problem04IdentificationFigures] = problem04Identification();
saveExperiment(fullfile(problemDirectory,'results'),'identification',problem04IdentificationResults,problem04IdentificationFigures);
[problem04MarkovResults, problem04MarkovFigures] = problem04Markov(problem04IdentificationResults);
saveExperiment(fullfile(problemDirectory,'results'),'markov',problem04MarkovResults,problem04MarkovFigures);
