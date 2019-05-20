function class_calc = knnclass(neighbors_class,neighbors_distance,num_class,K)

% find class on the basis of K neighbors
for g = 1:num_class
    freq(g) = length(find(neighbors_class == g));
end

unique_class = find(freq == max(freq));

if length(unique_class) == 1
    [M,class_calc] = max(freq);
else
    mean_dist = ones(num_class,1).*max(neighbors_distance);
    for g = 1:length(unique_class)
        in = find(neighbors_class == unique_class(g));
        mean_dist(unique_class(g)) = mean(neighbors_distance(in));
    end
    [m,class_calc] = min(mean_dist);
end