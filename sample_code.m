%%
%% load images and match files for the first example
%%

I1 = imread('library1.jpg');
I2 = imread('library2.jpg');
matches = load('library_matches.txt'); 
% this is a N x 4 file where the first two numbers of each row
% are coordinates of corners in the first image and the last two
% are coordinates of corresponding corners in the second image: 
% matches(i,1:2) is a point in the first image
% matches(i,3:4) is a corresponding point in the second image

N = size(matches,1);

%%

I1_two = rgb2gray(I1);
I1_two = im2double(I1_two);
I2_two = rgb2gray(I2);
I2_two = im2double(I2_two);

feat_points_left = harris(I1_two,1); %before it was 1
feat_points_left(:,1:100) = 0;
[c_left,r_left,cim_left] = find(feat_points_left > .04); %before it was .04
feat_points_right = harris(I2_two,1);
feat_points_right(:,300:end) = 0;
[c_right,r_right,cim_right] = find(feat_points_right > .04);

sift_points_left = find_sift(I1_two,[r_left,c_left,cim_left],2); %before it was 2
sift_points_right = find_sift(I2_two,[r_right,c_right,cim_right],2);

dist_matrix = dist2(sift_points_left,sift_points_right);

[matching_features_left,matching_features_right,distances] = find(dist_matrix < .03); %before it was .03

matches_sift = [r_left(matching_features_left(:)), c_left(matching_features_left(:)),r_right(matching_features_right(:)),c_right(matching_features_right(:))];








%%
%% display two images side-by-side with matches
%% this code is to help you visualize the matches, you don't need
%% to use it to produce the results for the assignment
%%
imshow([I1 I2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');

%%

imshow([I1 I2]); hold on;
plot(matches_sift(:,1), matches_sift(:,2), '+r');
plot(matches_sift(:,3)+size(I1,2), matches_sift(:,4), '+r');
line([matches_sift(:,1) matches_sift(:,3) + size(I1,2)]', matches_sift(:,[2 4])', 'Color', 'r');


%%
%% display second image with epipolar lines reprojected 
%% from the first image
%%

% first, fit fundamental matrix to the matches
type=1;
F = fit_fundamental(matches,type); % this is a function that you should write

L = (F * [matches(:,1:2) ones(N,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);

num_inliers = length(pt_line_dist(abs(pt_line_dist) < 3));
mean_distances = sum(pt_line_dist.^2)/N; %mean squared distance 
%for RANSAC method: 24.226 with 25 inliers (110)  on houses
%for unnormalized: 26.4013 with 31 inliers (121)  on houses
%for normalized: .07 with 168 inliers  on houses

%for RANSAC method: 6.32 with 250 inliers (110)  on library
%for unnormalized: .1792 with all inliers (309)  on library
%for normalized: .0603 with all inliers  on library



closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
clf;
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

%%
N2 = size(matches_sift,1);
type=3;
F = fit_fundamental(matches_sift,type); % this is a function that you should write

L = (F * [matches_sift(:,1:2) ones(N2,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches_sift(:,3:4) ones(N2,1)],2);
closest_pt = matches_sift(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
clf;
imshow(I2); hold on;
plot(matches_sift(:,3), matches_sift(:,4), '+r');
line([matches_sift(:,3) closest_pt(:,1)]', [matches_sift(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

%%
%camera center stuff

P1 = load('library1_camera.txt');
P2 = load('library2_camera.txt');

[U,S,V]=svd(P1); 
camera_1_center = V(:,end); %camera center

% imshow(I1); hold on;
% plot(camera_1_center(1:2), '+r');

[U,S,V]=svd(P2); 
camera_2_center = V(:,end); %camera center
% figure();
% imshow(I2); hold on;
% plot(camera_2_center(1:2), '+r');

%%
%triangulation

new_3d_points = zeros(N,3);
for i = 1:N
    left_pixel = [matches(i,1) matches(i,2) 1];
    right_pixel = [matches(i,3) matches(i,4) 1];

    
    A = [left_pixel(1)*P1(3,1) left_pixel(1)*P1(3,2) -P1(1,1) -P1(1,2);...
        left_pixel(2)*P1(3,1) left_pixel(2)*P1(3,2) -P1(2,1) -P1(2,2) ;...
        right_pixel(1)*P2(3,1) right_pixel(1)*P2(3,2) -P2(1,1) -P2(1,2) ;...
        right_pixel(2)*P2(3,1) right_pixel(2)*P2(3,2) -P2(2,1) -P2(2,2)];
    
    [U,S,V]=svd(A);
    X = V(:,end);
    X = X/X(4);
    new_3d_points(i,1:3) = X(1:3);

    
end

camera_1_center = camera_1_center/camera_1_center(4);
camera_2_center = camera_2_center/camera_2_center(4);

scatter3(new_3d_points(:,1),new_3d_points(:,2),new_3d_points(:,3));
hold on
scatter3(camera_1_center(1),camera_1_center(2),camera_1_center(3),'*');
hold on
scatter3(camera_2_center(1),camera_2_center(2),camera_2_center(3),'*');
axis equal


