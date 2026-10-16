% Problem 3: model analysis, level PI control, and independent loop comparison.
problemDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(problemDirectory,'..','common'),problemDirectory);
[problem03AnalysisResults, problem03AnalysisFigures] = problem03Analysis();
saveExperiment(fullfile(problemDirectory,'results'),'analysis',problem03AnalysisResults,problem03AnalysisFigures);
[problem03LevelResults, problem03LevelFigures] = problem03LevelControl();
saveExperiment(fullfile(problemDirectory,'results'),'level_control',problem03LevelResults,problem03LevelFigures);
[problem03ComparisonResults, problem03ComparisonFigures] = problem03Comparison();
saveExperiment(fullfile(problemDirectory,'results'),'comparison',problem03ComparisonResults,problem03ComparisonFigures);
