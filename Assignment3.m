
%assignment 3, part 1 stitching two images together

imname_left = 'uttower_left.JPG';
imname_right = 'uttower_right.JPG';

fullim_left = imread(imname_left);
color_left = fullim_left;
fullim_left = rgb2gray(fullim_left);
fullim_left = im2double(fullim_left);
fullim_right = imread(imname_right);
color_right = fullim_right;
fullim_right = rgb2gray(fullim_right);
fullim_right = im2double(fullim_right);
%figure, imshow(fullim_right);

feat_points_left = harris(fullim_left,1); %before it was 2
feat_points_left(:,1:512) = 0;
[c_left,r_left,cim_left] = find(feat_points_left > .15); %before it was .15
feat_points_right = harris(fullim_right,1);
feat_points_right(:,512:end) = 0;
[c_right,r_right,cim_right] = find(feat_points_right > .15);
%figure, imshow(feat_points_left);
%figure, imshow(feat_points_right);
%show_all_circles(fullim_left,r_left,c_left,cim_left);
%figure();
%show_all_circles(fullim_right,r_right,c_right,cim_right);

sift_points_left = find_sift(fullim_left,[r_left,c_left,cim_left],2); %before it was 1.5
sift_points_right = find_sift(fullim_right,[r_right,c_right,cim_right],2);

dist_matrix = dist2(sift_points_left,sift_points_right);

[matching_features_left,matching_features_right,distances] = find(dist_matrix < .04); %previously .017
%%RANSAC portion
%%
points_to_fit = [r_left(matching_features_left(:)), c_left(matching_features_left(:)),r_right(matching_features_right(:)),c_right(matching_features_right(:))];
num_inliers = 0;
while true
    
    rand_nums = randsample(196,4);
    rand_inliers = zeros(4,4); 
    rand_inliers(1,:) = [r_left(matching_features_left(rand_nums(1))), c_left(matching_features_left(rand_nums(1))), r_right(matching_features_right(rand_nums(1))), c_right(matching_features_right(rand_nums(1)))];
    rand_inliers(2,:) = [r_left(matching_features_left(rand_nums(2))), c_left(matching_features_left(rand_nums(2))), r_right(matching_features_right(rand_nums(2))), c_right(matching_features_right(rand_nums(2)))];
    rand_inliers(3,:) = [r_left(matching_features_left(rand_nums(3))), c_left(matching_features_left(rand_nums(3))), r_right(matching_features_right(rand_nums(3))), c_right(matching_features_right(rand_nums(3)))];
    rand_inliers(4,:) = [r_left(matching_features_left(rand_nums(4))), c_left(matching_features_left(rand_nums(4))), r_right(matching_features_right(rand_nums(4))), c_right(matching_features_right(rand_nums(4)))];
    %rand_inliers format 4x[row_left,column_left,row_right,column_right]
    A = zeros(8,9);
    A(1,:) = [0,0,0,rand_inliers(1,1),rand_inliers(1,2),1,-(rand_inliers(1,4)*rand_inliers(1,1)),-(rand_inliers(1,4)*rand_inliers(1,2)),-(rand_inliers(1,4)*1)];
    A(2,:) = [rand_inliers(1,1),rand_inliers(1,2),1,0,0,0,-(rand_inliers(1,3)*rand_inliers(1,1)),-(rand_inliers(1,3)*rand_inliers(1,2)),-(rand_inliers(1,3)*1)];
    A(3,:) = [0,0,0,rand_inliers(2,1),rand_inliers(2,2),1,-(rand_inliers(2,4)*rand_inliers(2,1)),-(rand_inliers(2,4)*rand_inliers(2,2)),-(rand_inliers(2,4)*1)];
    A(4,:) = [rand_inliers(2,1),rand_inliers(2,2),1,0,0,0,-(rand_inliers(2,3)*rand_inliers(2,1)),-(rand_inliers(2,3)*rand_inliers(2,2)),-(rand_inliers(2,3)*1)];
    A(5,:) = [0,0,0,rand_inliers(3,1),rand_inliers(3,2),1,-(rand_inliers(3,4)*rand_inliers(3,1)),-(rand_inliers(3,4)*rand_inliers(3,2)),-(rand_inliers(3,4)*1)];
    A(6,:) = [rand_inliers(3,1),rand_inliers(3,2),1,0,0,0,-(rand_inliers(3,3)*rand_inliers(3,1)),-(rand_inliers(3,3)*rand_inliers(3,2)),-(rand_inliers(3,3)*1)];
    A(7,:) = [0,0,0,rand_inliers(4,1),rand_inliers(4,2),1,-(rand_inliers(4,4)*rand_inliers(4,1)),-(rand_inliers(4,4)*rand_inliers(4,2)),-(rand_inliers(4,4)*1)];
    A(8,:) = [rand_inliers(4,1),rand_inliers(4,2),1,0,0,0,-(rand_inliers(4,3)*rand_inliers(4,1)),-(rand_inliers(4,3)*rand_inliers(4,2)),-(rand_inliers(4,3)*1)];
    %solve Ah=0 for h using least squares which in this case is AX = 0
    [U,S,V]=svd(A); 
    X = V(:,end); %transformation vector of h1-h9
    H_mat = zeros(3,3);
    H_mat(1,:) = X(1:3);
    H_mat(2,:) = X(4:6);
    H_mat(3,:) = X(7:9);
    %find l2 norm of (p2 - transform p1) 
    ssd_output = zeros(196,1);
    inliers_to_display = zeros(196,4);
    for i = 1:196
        transform_p1 = H_mat*transpose([points_to_fit(i,1:2),1]);
        transform_p1 = transform_p1./transform_p1(3);
        curr_point_dist = points_to_fit(i,3:4) - transpose(transform_p1(1:2));
        ssd = sum(curr_point_dist(:).^2);
        ssd_output(i) = ssd;
        if ssd < .002
            inliers_to_display(i,:) = points_to_fit(i,:);
            num_inliers = num_inliers+1;
        end
    end
    if num_inliers > 15
        break;
    else
        num_inliers = 0;
    end
end

%%
%number of inliers is about 21 out of 196 matching features
%residual avergae is .002 for those inliers



inliers_to_display = inliers_to_display(any(inliers_to_display,2),:);
imshow([fullim_left fullim_right]); hold on;
plot(inliers_to_display(:,1), inliers_to_display(:,2), '+r');
plot(inliers_to_display(:,3)+size(fullim_left,2), inliers_to_display(:,4), '+r');
line([inliers_to_display(:,1) inliers_to_display(:,3) + size(fullim_left,2)]', inliers_to_display(:,[2 4])', 'Color', 'r');

%%
   

    [h_inliers w_inliers] = size(inliers_to_display);
    rand_nums_2 = randsample(h_inliers,4);
    %T = maketform('projective', [inliers_to_display(rand_nums_2,1) inliers_to_display(rand_nums_2,2)], [inliers_to_display(rand_nums_2,3) inliers_to_display(rand_nums_2,4)]);
    T = maketform('projective',[rand_inliers(:,1) rand_inliers(:,2)],[rand_inliers(:,3) rand_inliers(:,4)]);
    %new_img_2 = imtransform(fullim_left,T);
    tform = projective2d(transpose(H_mat));
    new_img = imwarp(fullim_left,tform);


    [im2t,xdataim2t,ydataim2t]=imtransform(fullim_left,T);
    xdataout=[min(1,xdataim2t(1)) max(size(fullim_right,2),xdataim2t(2))];
    ydataout=[min(1,ydataim2t(1)) max(size(fullim_right,1),ydataim2t(2))];
    im2t=imtransform(fullim_left,T,'XData',xdataout,'YData',ydataout);
    im1t=imtransform(fullim_right,maketform('affine',eye(3)),'XData',xdataout,'YData',ydataout);
    ims=im1t/2+im2t/2;
    figure, imshow(ims);

    %T2 = maketform('projective', [inliers_to_display(rand_nums_2,1) inliers_to_display(rand_nums_2,2)], [inliers_to_display(rand_nums_2,3) inliers_to_display(rand_nums_2,4)]);
    T2 = maketform('projective',[rand_inliers(:,1) rand_inliers(:,2)],[rand_inliers(:,3) rand_inliers(:,4)]);
    [im2t,xdataim2t,ydataim2t]=imtransform(color_left,T2);
    xdataout=[min(1,xdataim2t(1)) max(size(color_right,2),xdataim2t(2))];
    ydataout=[min(1,ydataim2t(1)) max(size(color_right,1),ydataim2t(2))];
    im2t=imtransform(color_left,T2,'XData',xdataout,'YData',ydataout);
    im1t=imtransform(color_right,maketform('affine',eye(3)),'XData',xdataout,'YData',ydataout);
    ims=im1t/2+im2t/2;
    figure, imshow(ims);

%RANSAC loop: 
% 1. Select four feature pairs (at random) 
% 2. Compute homography H (exact) 
% 3.Compute inliers where  SSD(pi’, Hpi)< ε
% 4.Keep largest set of inliers 
% 5.Re-compute least-squares H estimate on all of the inliers



