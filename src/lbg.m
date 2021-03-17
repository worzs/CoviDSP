function new_centroids = lbg(samples, M_max, step_size, error_threshold)
    [sample_num, feature_num] = size(samples);    % Time & MFCCs
    codebook_index = zeros(sample_num, 1); 
    
    old_centroids = reshape(mean(samples, 1), 1, feature_num); %initial centroid = mean of samples(training data)
    
    M = 1;
    old_error = sum(vecnorm(samples - old_centroids, 2, 2)); % compute the first error
    fprintf("M = %d, error = %f \n", M, old_error);
    
    while(M < M_max)                              % required M clusters

        new_centroids = zeros(2 * M, feature_num);% initialize centroid array
                                                  % splitting centroids
        for i = 1: M
            new_centroids(2 * i - 1, :) = old_centroids(i, :) * (1 + step_size);
            new_centroids(2 * i, :) = old_centroids(i, :) * (1 - step_size);
        end

        M = M * 2;
        while (true)

            cnts = zeros(M, 1);                       % # of training vector assigned to each clusters
            for i = 1: sample_num
                distance = zeros(M, 1);
                for j = 1: M
                    distance(j) = norm(samples(i, :) - new_centroids(j, :), 2); % Euclidean distance between MFCC and centroid
                end
                [~, ind] = min(distance);
                codebook_index(i) = ind(1);           % cluster index of the training vector
                cnts(ind(1)) = cnts(ind(1)) + 1;
            end 

            old_centroids = zeros(M, feature_num);    
            for i = 1: sample_num
                old_centroids(codebook_index(i), :) = old_centroids(codebook_index(i), :) + samples(i, :);
            end

            for i = 1: M
                if (cnts(i) ~= 0)
                    old_centroids(i, :) = old_centroids(i, :) / cnts(i);
                end
            end

            new_centroids = old_centroids;

            new_error = 0;
            for i = 1: sample_num
                new_error = new_error + norm(samples(i, :) - old_centroids(codebook_index(i), :), 2);
            end

            fprintf("M = %d, error = %f \n", M, new_error);
            if (new_error == 0)
                break
            else
                relative_error = abs(new_error - old_error) / old_error;
                old_error = new_error;
                if (relative_error < error_threshold)
                    break;
                end
            end
        end
    end

    %new_centroids(all(~new_centroids, 2), :) = []; % Remove zero rows