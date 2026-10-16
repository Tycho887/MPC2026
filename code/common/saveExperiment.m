function saveExperiment(directory, name, result, figures)
%SAVEEXPERIMENT Save explicit data and figure handles as report artifacts.
% Existing artifacts with the same experiment name are replaced on rerun.
if ~isfolder(directory), mkdir(directory); end
save(fullfile(directory,[name '.mat']),'result');
for k = 1:numel(figures)
    stem = fullfile(directory,sprintf('%s_%02d',name,k));
    exportgraphics(figures(k),[stem '.png'],'Resolution',200);
    exportgraphics(figures(k),[stem '.pdf'],'ContentType','vector');
end
end
