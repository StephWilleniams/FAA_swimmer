% Get the jump arrays.

fi_counter = 0;

for fi_i = [1,40]
    for fi_j = [1,12]
        for k = 1:10

            fi_counter = fi_counter +1;
            pass = load(['outputs_extrema/' num2str(fi_j) '_' num2str(fi_i) '_' num2str(fi_counter) '/output_passive.txt']); % Load in the data

            Nparts = 40; % Number of particles being checked.
            TperP = length(pass(:,1))/Nparts; % Timesteps per particle.
            passN = zeros(TperP,5); % Store the the nth particle's data.
            start_end_id = zeros(TperP,1); % store for the start/end data (based on if F_x is non-zero).
            jump_count = 0; % Counter for the number of jumps.

            % Count the jumps to preallocate an array for them.
            for n = 1:40
                a = find(pass(:,2) == n); % Get the indicies of particle n.
                passN = pass(a,:); % Get the trajectory of particle n.
                start_end_id(1:end-1) = passN(1:end-1,5)-passN(2:end,5); % If -1 jump start, if 1 jump end.
                starts = find(start_end_id == -1); % Get the jump starts.
                ends = find(start_end_id == 1); % Get the jump ends.
                mi = min(length(starts),length(ends)); % Get the minimum of starts/ends.
                starts = starts(1:mi); % Get the minimum of starts.
                ends = ends(1:mi); % Get the minimum of ends.
                jump_count = jump_count + length(starts); % Count the number of jumps for this particle.
            end

            jumps = zeros(jump_count,6); % Preallocate space for all the jumps calculated.
            counter = 0; % Jump number counter.

            % Fill the preallocated array.
            for n = 1:40
                % As in previous section
                a = find(pass(:,2) == n);
                passN = pass(a,:);
                start_end_id(1:end-1) = passN(1:end-1,5)-passN(2:end,5); % If -1 jump start, if 1 jump end.
                starts = find(start_end_id == -1); % Get the jump starts.
                jumps = jumps + length(starts);
                ends = find(start_end_id == 1); % Get the jump ends.
                mi = min(length(starts),length(ends));
                starts = starts(1:mi);
                ends = ends(1:mi);

                % Extract the info on the jumps
                for i = 1:length(starts)
                    counter = counter+1;
                    jumps(counter,1) = passN(starts(i),1);
                    jumps(counter,2) = passN(ends(i),1);
                    jumps(counter,3) = passN(starts(i),3);
                    jumps(counter,4) = passN(ends(i),3);
                    jumps(counter,5) = passN(starts(i),4);
                    jumps(counter,6) = passN(ends(i),4);
                end
            end

            save(['outputs_extrema/jumps/' num2str(fi_j) '_' num2str(fi_i) '_' num2str(fi_counter) '_jumps.mat'],"jumps");

        end
    end
end