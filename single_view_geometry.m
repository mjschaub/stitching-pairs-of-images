%%part 3 code
%%

imname = 'north_quad.jpg';

fullim = imread(imname);
imshow(fullim);

%%
vp_x = getVanishingPoint_shell(fullim);
%%
vp_y = getVanishingPoint_shell(fullim);
%%
vp_z = getVanishingPoint_shell(fullim);


%%
original_vp_x = [-2.279016182971489e+02;2.072417796449178e+02;1]; %vp_x
original_vp_y = [1.411587433926172e+03;2.359969337938032e+02;1];
original_vp_z = [5.914377363320575e+02;6.987089155443199e+03;1];

%%

vp_y = original_vp_y;
vp_z = original_vp_z;
vp_x = original_vp_x;


%%
vp_y = original_vp_y;
vp_z = original_vp_z;
vp_x = original_vp_x;
figure();
imagesc(fullim);
hold on;
plot([vp_x(1) vp_y(1)], [vp_x(2) vp_y(2)]);
axis image;

horizon_line = real(cross(vp_x', vp_y'));

length = sqrt(horizon_line(1)^2 + horizon_line(2)^2);
horizon_line = horizon_line/length;

%eq = horizon_line(1)*x + horizon_line(2)*y + horizon_line(3) = 0
%horizon line: -.0175x+.9998y-211.2065 = 0

%z -> x
%x -> y
%y -> z
%%
%optical center u and v
%focal length f
% set left most vanishing point as z, right most as x, and vertical as y
vp_y = original_vp_z;
vp_z = original_vp_x;
vp_x = original_vp_y;

temp_1_x = vp_x(1)+vp_y(1);
temp_2_x = vp_y(1)+vp_z(1);
temp_3_x = vp_z(1)+vp_x(1);
temp_1_y = vp_x(2)+vp_y(2);
temp_2_y = vp_y(2)+vp_z(2);
temp_3_y = vp_z(2)+vp_x(2);
temp_xy = vp_x(1)*vp_y(1)+vp_x(2)*vp_y(2);
temp_yz = vp_z(1)*vp_y(1)+vp_z(2)*vp_y(2);
temp_xz = vp_x(1)*vp_z(1)+vp_x(2)*vp_z(2);

%solve for image center ax+by+c=0
syms u_0 v_0;
[final_u, final_v] = solve(-(temp_1_x)*u_0 - (temp_1_y)*v_0 + temp_xy == -(temp_2_x)*u_0 -(temp_2_y)*v_0 + temp_yz,...
-(temp_3_x)*u_0 - (temp_3_y)*v_0 + temp_xz == -(temp_2_x)*u_0 - (temp_2_y)*v_0 + temp_yz);

u = double(final_u);
v = double(final_v);
%focal length solving
eq_f = (u - vp_x(1))*(u - vp_y(1)) + (v - vp_x(2))*(v - vp_y(2));
f = sqrt(-eq_f);

K = [f      0       u;
     0      f       v;
     0      0       1];    
 
 
 %focal_length = 805.4438
 %image center = (708.3552,320.982)
 %%
%rotation matrix  vp_i = K*R 

R_two = K \ [vp_x(1) vp_x(2) vp_x(3);vp_y(1) vp_y(2) vp_y(3);vp_z(1) vp_z(2) vp_z(3)];
r_x = K\vp_x;
r_y = K\vp_y;
r_z = K\vp_z;

unit_rx = sqrt(sumsqr(r_x)); 
unit_ry = sqrt(sumsqr(r_y)); 
unit_rz = sqrt(sumsqr(r_z)); 

r_x = r_x/unit_rx;
r_y = r_y/unit_ry;
r_z = r_z/unit_rz;

R = [r_x r_y r_z];

%R = [0.655625547050926,-0.0174098133942597,-0.754885448562728;-0.0792318444389579,0.992628949695479,-0.0917065049669339;0.750917743085736,0.119936093916829,0.649413486536120]



%%
%%measuring height of objects

bottom_guy = [627 515 1]; %bottom pixel of guy
top_guy = [627 467 1]; %top pixel of guy
height_of_person = 66; %5ft 6 in in inches is 66 or 72 for 6 feet

%%
%CSL Building
bottom_CSL = [600 328 1];
top_CSL = [600 55 1];

line_between = real(cross(bottom_guy', bottom_CSL'));
vanishing_line = real(cross(line_between, horizon_line));
vanishing_line = vanishing_line/vanishing_line(3);

line_temp = real(cross(vanishing_line', top_guy'));
height_line = real(cross(top_CSL', bottom_CSL'));
t = real(cross(line_temp', height_line'));
t = t/t(3);

building_height = height_of_person * sqrt(sumsqr(vp_z' - t)) * sqrt(sumsqr(top_CSL - bottom_CSL)) / sqrt(sumsqr(t - bottom_CSL)) / sqrt(sumsqr(vp_z' - top_CSL));
%got 1025.3 inches or 85.4 feet for 5ft 6in
%got 1118.5 inches or 93.2 feet for 6ft

figure();
imagesc(fullim);
hold on; 
%plot([vanishing_line(1) bottom_guy(1)], [vanishing_line(2) bottom_guy(2)], 'r');
plot([top_guy(1) bottom_guy(1)], [top_guy(2) bottom_guy(2)], 'r');
%plot([vanishing_line(1) top_guy(1)], [vanishing_line(2) top_guy(2)], 'r');
%plot([vanishing_line(1) t(1)], [vanishing_line(2) t(2)], 'y');
%plot([bottom_CSL(1) t(1)], [bottom_CSL(2) t(2)], 'y');
%plot([bottom_CSL(1) vanishing_line(1)], [bottom_CSL(2) vanishing_line(2)], 'y');
plot([bottom_CSL(1) top_CSL(1)], [bottom_CSL(2) top_CSL(2)], 'g');
axis equal;
axis image;

%%
%spike statue
bottom_spike = [601 475 1];
top_spike = [601 254 1];

line_between = real(cross(bottom_guy', bottom_spike'));
vanishing_line = real(cross(line_between, horizon_line));
vanishing_line = vanishing_line/vanishing_line(3);

line_temp = real(cross(vanishing_line', top_guy'));
height_line = real(cross(top_spike', bottom_spike'));
t = real(cross(line_temp', height_line'));
t = t/t(3);

spike_height = height_of_person * sqrt(sumsqr(vp_z' - t)) * sqrt(sumsqr(top_spike - bottom_spike)) / sqrt(sumsqr(t - bottom_spike)) / sqrt(sumsqr(vp_z' - top_spike));
%spike height is 363.6 inches or 30.3 feet for 5ft 6in
%spike height is 396.66 inches or 33 feet for 6ft

figure();
imagesc(fullim);
hold on; 
plot([top_guy(1) bottom_guy(1)], [top_guy(2) bottom_guy(2)], 'r');
plot([bottom_spike(1) top_spike(1)], [bottom_spike(2) top_spike(2)], 'g');
axis equal;
axis image;

%%
%lamp posts
top_lamp = [1000 368 1];
bottom_lamp = [1000 483 1];

line_between = real(cross(bottom_guy', bottom_lamp'));
vanishing_line = real(cross(line_between, horizon_line));
vanishing_line = vanishing_line/vanishing_line(3);

line_temp = real(cross(vanishing_line', top_guy'));
height_line = real(cross(top_lamp', bottom_lamp'));
t = real(cross(line_temp', height_line'));
t = t/t(3);

lamp_height = height_of_person * sqrt(sumsqr(vp_z' - t)) * sqrt(sumsqr(top_lamp - bottom_lamp)) / sqrt(sumsqr(t - bottom_lamp)) / sqrt(sumsqr(vp_z' - top_lamp));
%lamp height is 183.8 inches or 15.3 feet for 5ft 6in
%lamp height is 200.51 inches or 16.7 feet for 6ft


figure();
imagesc(fullim);
hold on; 
plot([top_guy(1) bottom_guy(1)], [top_guy(2) bottom_guy(2)], 'r');
plot([bottom_lamp(1) top_lamp(1)], [bottom_lamp(2) top_lamp(2)], 'g');
axis equal;
axis image;

