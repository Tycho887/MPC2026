function [result, figures] = problem04Markov(identified)
% Migrated experiment; shared physical models live in common/.
figures = gobjects(0);
result = markovParameters(identified.G,4,60);
result.identification = identified;
H = result.H; t_k = result.t; Ts = result.Ts;
p_out = size(H,1); m_in = size(H,2);
%% 4. Plot Markov Parameters (Impulse Responses)
figures(end+1) = experimentFigure('Name', 'Problem 4.6: Discrete-Time Markov Parameters', ...
    'Color', 'w', 'Position', [150, 150, 950, 650]);

input_names  = {'F_1', 'F_2'};
output_names = {'h_1', 'h_2'};

for i_out = 1:p_out
    for j_in = 1:m_in
        subplot(p_out, m_in, (i_out - 1) * m_in + j_in);
        stem(t_k, squeeze(H(i_out, j_in, :)), 'filled', ...
            'MarkerSize', 4, 'LineWidth', 1.2, 'Color', [0.1, 0.4, 0.8]);
        grid on;
        xlabel('Time (s)');
        ylabel(sprintf('H_{%d%d}(k)', i_out, j_in));
        title(sprintf('%s \\rightarrow %s (T_s = %.1f s)', ...
            input_names{j_in}, output_names{i_out}, Ts));
    end
end



end
