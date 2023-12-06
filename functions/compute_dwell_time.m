function [dwell_times] = compute_dwell_time(kt, plot_logic)
markov_chain = kt; % Sequence of states

% Initialize variables
unique_states = unique(markov_chain);  % Get unique states
dwell_times = cell(length(unique_states), 1);  % Initialize dwell times as a cell array

% Compute dwell times
current_state = markov_chain(1);  % Initialize with the first state
dwell_time = 1;

for i = 2:length(markov_chain)
    if markov_chain(i) == current_state
        % If the current state is the same as the previous state, increment dwell time
        dwell_time = dwell_time + 1;
    else
        % If the current state is different, store the dwell time and reset counter
        state_idx = find(unique_states == current_state);  % Find the index of the current state
        dwell_times{state_idx} = [dwell_times{state_idx}, dwell_time];
        dwell_time = 1;
        current_state = markov_chain(i);  % Update the current state
    end
end

% Add the last dwell time to the list
state_idx = find(unique_states == current_state);  % Find the index of the last state
dwell_times{state_idx} = [dwell_times{state_idx}, dwell_time];

if isempty(plot_logic)~=1
% Plot the dwell time distribution for each state
for i = 1:length(unique_states)
    state = unique_states(i);
    dwell_time_data = dwell_times{i};
    
    figure;
    histogram(dwell_time_data, 'Normalization', 'probability');
    xlabel('Dwell Time');
    ylabel('Probability');
    title(['Dwell Time Distribution for State ' num2str(state)]);
end
end

end
