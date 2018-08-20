function fit_fund = fit_fundamental(matches,type)
%%
%% Fit fundamental matrix to the matches
%% this is a N x 4 file where the first two numbers of each row
%% are coordinates of corners in the first image and the last two
%% are coordinates of corresponding corners in the second image: 
%% matches(i,1:2) is a point in the first image
%% matches(i,3:4) is a corresponding point in the second image
%%
fprintf('Running fit_fundamental\n');

N = size(matches,1);

if type == 0 %unnormalized

    %unnormalized
    rand_nums = randsample(N,8);

    A = zeros(8,9);
    ones_vec = ones(8,1);
    A(:,:) = [matches(rand_nums(:),3).*matches(rand_nums(:),1),matches(rand_nums(:),3).*matches(rand_nums(:),2),matches(rand_nums(:),3),matches(rand_nums(:),4).*matches(rand_nums(:),1),matches(rand_nums(:),4).*matches(rand_nums(:),2),matches(rand_nums(:),4),matches(rand_nums(:),1),matches(rand_nums(:),2),ones_vec];
    
    
    temp_A = zeros(N,9);
    one_vec = ones(N,1);
    temp_A(:,:) = [matches(:,3).*matches(:,1),matches(:,3).*matches(:,2),matches(:,3),matches(:,4).*matches(:,1),matches(:,4).*matches(:,2),matches(:,4),matches(:,1),matches(:,2),one_vec];

    [U,S,V]=svd(temp_A); 
    X = V(:,end); %transformation vector of f1-f9
    F = reshape(X,[3 3])';
    [U_2,S_2,V_2] = svd(F);
    S_2(3,3) = 0;
    fit_fund = U_2*S_2*V_2';
    

elseif type == 1

    % how to normalize
    matches_mean = mean(matches,1);
    mean_x = ones(N,1)*matches_mean(1);
    mean_y = ones(N,1)*matches_mean(2);
    mean_x_2 = ones(N,1)*matches_mean(3);
    mean_y_2 = ones(N,1)*matches_mean(4);
    mean_first_img = zeros(N,2);
    
    mean_first_img(:,:) = [mean_x mean_y];
    dists_1 = sqrt(sum((matches(:,1:2) - mean_first_img).^2,2));
    mean_dist_1 = mean(dists_1);
    
    mean_second_img(:,:) = [mean_x_2 mean_y_2];
    dists_2 = sqrt(sum((matches(:,3:4) - mean_second_img).^2,2));
    mean_dist_2 = mean(dists_2);
    
    T = [sqrt(2)/mean_dist_1 0 -sqrt(2)/mean_dist_1*matches_mean(1); 0 sqrt(2)/mean_dist_1 -sqrt(2)/mean_dist_1*matches_mean(2); 0 0 1];  %points 1 and 2
    T_2 = [sqrt(2)/mean_dist_2 0 -sqrt(2)/mean_dist_2*matches_mean(3); 0 sqrt(2)/mean_dist_2 -sqrt(2)/mean_dist_2*matches_mean(4); 0 0 1]; %points 3 and 4
    new_matches = zeros(N,4);
    for i = 1:N
        
        normalized_left_coord = T*transpose([matches(i,1:2),1]);
        normalized_right_coord = T_2*transpose([matches(i,3:4),1]);
        normalized_left_coord = normalized_left_coord./normalized_left_coord(3);
        normalized_right_coord = normalized_right_coord./normalized_right_coord(3);
        new_matches(i,:) = [normalized_left_coord(1),normalized_left_coord(2),normalized_right_coord(1),normalized_right_coord(2)];
        
    end
    
    A_normalized = zeros(N,9); %still have to do, then get F and change it back
    one_vec = ones(N,1);
    A_normalized(:,:) = [new_matches(:,3).*new_matches(:,1),new_matches(:,3).*new_matches(:,2),new_matches(:,3),new_matches(:,4).*new_matches(:,1),new_matches(:,4).*new_matches(:,2),new_matches(:,4),new_matches(:,1),new_matches(:,2),one_vec];
    
    [U,S,V]=svd(A_normalized); 
    X_2 = V(:,end); %transformation vector of f1-f9
    F_2 = reshape(X_2,[3 3])';
    [U_2,S_2,V_2] = svd(F_2);
    S_2(3,3) = 0;
    F_temp = U_2*S_2*V_2';
    
    fit_fund = T_2'* F_temp * T;
    
    
else
    
    matches_mean = mean(matches,1);
    mean_x = ones(N,1)*matches_mean(1);
    mean_y = ones(N,1)*matches_mean(2);
    mean_x_2 = ones(N,1)*matches_mean(3);
    mean_y_2 = ones(N,1)*matches_mean(4);
    mean_first_img = zeros(N,2);
    
    mean_first_img(:,:) = [mean_x mean_y];
    dists_1 = sqrt(sum((matches(:,1:2) - mean_first_img).^2,2));
    mean_dist_1 = mean(dists_1);
    
    mean_second_img(:,:) = [mean_x_2 mean_y_2];
    dists_2 = sqrt(sum((matches(:,3:4) - mean_second_img).^2,2));
    mean_dist_2 = mean(dists_2);
    
    T = [sqrt(2)/mean_dist_1 0 -sqrt(2)/mean_dist_1*matches_mean(1); 0 sqrt(2)/mean_dist_1 -sqrt(2)/mean_dist_1*matches_mean(2); 0 0 1];  %points 1 and 2
    T_2 = [sqrt(2)/mean_dist_2 0 -sqrt(2)/mean_dist_2*matches_mean(3); 0 sqrt(2)/mean_dist_2 -sqrt(2)/mean_dist_2*matches_mean(4); 0 0 1]; %points 3 and 4
    new_matches = zeros(N,4);
    for i = 1:N
        
        normalized_left_coord = T*transpose([matches(i,1:2),1]);
        normalized_right_coord = T_2*transpose([matches(i,3:4),1]);
        normalized_left_coord = normalized_left_coord./normalized_left_coord(3);
        normalized_right_coord = normalized_right_coord./normalized_right_coord(3);
        new_matches(i,:) = [normalized_left_coord(1),normalized_left_coord(2),normalized_right_coord(1),normalized_right_coord(2)];
        
    end
    
    %A_normalized = zeros(N,9); %still have to do, then get F and change it back
    %one_vec = ones(N,1);
    one_vec_rand = ones(8,1);
    rand_nums = randsample(N,8);
    A_random = zeros(8,9);
    A_random(:,:) = [new_matches(rand_nums,3).*new_matches(rand_nums,1),new_matches(rand_nums,3).*new_matches(rand_nums,2),new_matches(rand_nums,3),new_matches(rand_nums,4).*new_matches(rand_nums,1),new_matches(rand_nums,4).*new_matches(rand_nums,2),new_matches(rand_nums,4),new_matches(rand_nums,1),new_matches(rand_nums,2),one_vec_rand];
    %A_normalized(:,:) = [new_matches(:,3).*new_matches(:,1),new_matches(:,3).*new_matches(:,2),new_matches(:,3),new_matches(:,4).*new_matches(:,1),new_matches(:,4).*new_matches(:,2),new_matches(:,4),new_matches(:,1),new_matches(:,2),one_vec];
    
    [U,S,V]=svd(A_random); %A_normalized
    X_2 = V(:,end); %transformation vector of f1-f9
    F_2 = reshape(X_2,[3 3])';
    [U_2,S_2,V_2] = svd(F_2);
    S_2(3,3) = 0;
    F_temp = U_2*S_2*V_2';
    
    
    
    fit_fund = T_2'* F_temp * T;
    
    
    
    
    
end


