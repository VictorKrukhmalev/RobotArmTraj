% the script on recovering kinematics of 2-linked manipulator kinematic
% motion on the basis of sensors data

clc
clear
close all
% open srd1.xls accelerometers on link1 data, preliminary saved at .mat

accs_data = open('srd1_accs.mat');
% accs_data = open('srd2_accs.mat');

% link 1 accelerations
a1_l1 = accs_data.l1a1;
a2_l1 = accs_data.l1a2;
a3_l1 = accs_data.l1a3;
a4_l1 = accs_data.l1a4;%this is IMU accelerations

% link 1 accelerations
a1_l2 = accs_data.l2a1;
a2_l2 = accs_data.l2a2;
a3_l2 = accs_data.l2a3;
a4_l2 = accs_data.l2a4;%this is IMU accelerations

%initial data on task
dt=0.01;%time step, s, considering rate 100Hz
a_gravity_gl = [0,0,-9.8];%3d vector of gravity

% link 1 data
link1_initpose = [178,-38.96,-144.15];%initial euler angles of link 1 in ground frame
l1_len = 0.4;%length of link1 body
l1_senslen = 0.2;%length where sensors located on link 1
V_l1_gl = [0;0;0];%initial velocity of link 1 is zero
V_l1 = [0;0;0];%initial velocity of link 1 is zero
%assuming that ZERO pose mean than both link1 and link2 are on X axis
l1_senspoint = [0.2, 0., 0.];
l1_endpoint = [0.4, 0., 0.];

%link 2 init dat
link2_initpose = [-90.50 63.73 -99.97];
l2_len = 0.4;
l2_senslen = 0.2;
V_l2_gl = [0;0;0];%initial velocity of link 2 in base frame is zero
V_l2_onl1 = [0;0;0];%initial velocity of link 2 in link1 reference frame is zero
l2_senspoint_zero = [0.2, 0., 0.];%zero pose of l2 sensor point on link1
l2_endpoint_zero = [0.4, 0., 0.];%zero pose of l2 end point on link1

% all rotations we shall do with quaternions, so we avoid 0-360 issues

%first we convert initial pose to quaternion
%as far as matlab function do it ZYX, we flip init_pose vector from right
%to left, converting from degrees to radians

% l1_quat_total = eul2quat(deg2rad(fliplr(link1_initpose)),'ZYX')
link1_initpose_rads = deg2rad(link1_initpose);
l1_quat_total = angle2quat(link1_initpose_rads(1),link1_initpose_rads(2),link1_initpose_rads(3),'XYZ');

%according to task we have initial pose in ground frame
link2_initpose_rads = deg2rad(link2_initpose);
l2_quat_total_base = angle2quat(link2_initpose_rads(1),link2_initpose_rads(2),link2_initpose_rads(3),'XYZ');
%we find angles between links continuing XYZ sequence paradigm
%rotation of link 2 at link1 (in reference to link 1 frame)
l2_quat_onl1_total = quatmultiply(l2_quat_total_base, quatinv(l1_quat_total))


% and here we go on our main loop
angl1_xlog = [];
angl1_ylog = [];
angl1_zlog = [];
omegal1_log = [ ];%this is for link 1
epsilon_log = [ ];%this is for link 1

angl2onl1_xlog = []; 
angl2onl1_ylog = []; 
angl2onl1_zlog = []; 

omega_l2onl1_log = [];
epsilon_l2onl1_log = [];

l2_base_angles_log = [];

for step = 1:length(a1_l1)
% for step = 1:2
    display(step)
    l1_quat_total;
    %
%     euler_l1 = quat2angle(l1_quat_total,'XYZ')
    [angl1_x, angl1_y, angl1_z]  = quat2angle(l1_quat_total,'XYZ');
    angl1_xlog = [angl1_xlog angl1_x]; 
    angl1_ylog = [angl1_ylog angl1_y]; 
    angl1_zlog = [angl1_zlog angl1_z]; 
    
    l1_senspoint_gl = quatrotate(l1_quat_total, l1_senspoint);
    l1_endpoint_gl = quatrotate(l1_quat_total, l1_endpoint);
    
    % as far as 4 accelerometers are located in circle we have to consider
    % 90degrees rotation of every sensors around the X axis
    %  and we neglect 2 sm shift of sensor from
    % link axis - in the result of bringin 4 acceleration to one point
    % resulting error of this neglection should be eliminated    
    a1onlink_l1 = a1_l1(step,:)';%if look along Xaxis - this is top sensor
    a2onlink_l1 = rotx(-90)*a2_l1(step,:)';%if look along Xaxis - this is right sensor
    a3onlink_l1 = rotx(-270)*a3_l1(step,:)';%left sensor
    a4onlink_l1 = rotx(-180)*a4_l1(step,:)';%bottom sensor
        
    % we have to rotate accs from local frame to global frame,
    % this is direct rotation: from link1 frame to global
    % frame
    a1onlink_l1_gl = quatrotate((l1_quat_total), a1onlink_l1');
    a2onlink_l1_gl = quatrotate((l1_quat_total), a2onlink_l1');
    a3onlink_l1_gl = quatrotate((l1_quat_total), a3onlink_l1');
    a4onlink_l1_gl = quatrotate((l1_quat_total), a4onlink_l1');
    
    %now eliminate gravity out of sensors data(the frame is global)
    a1motion_l1_gl = a1onlink_l1_gl - a_gravity_gl;
    a2motion_l1_gl = a2onlink_l1_gl - a_gravity_gl;
    a3motion_l1_gl = a3onlink_l1_gl - a_gravity_gl;
    a4motion_l1_gl = a4onlink_l1_gl - a_gravity_gl;
    %and get pseudo-equilibrium acceleration
    %TODO: dont we gonna do it on L1 frame?
    a_l1_gl = [0;0;0];    
    a_l1_gl(1) = norm([a1motion_l1_gl(1);a2motion_l1_gl(1);a3motion_l1_gl(1);a4motion_l1_gl(1) ]);
    a_l1_gl(2) = norm([a1motion_l1_gl(2);a2motion_l1_gl(2);a3motion_l1_gl(2);a4motion_l1_gl(2) ]);
    a_l1_gl(3) = norm([a1motion_l1_gl(3);a2motion_l1_gl(3);a3motion_l1_gl(3);a4motion_l1_gl(3) ]);

    dV_l1_gl = a_l1_gl*dt;
    V_l1_gl = V_l1_gl + dV_l1_gl;
    % go to local frame (because quaternion 'delta-summation' is only for
    %local) - from global we go with inverse quaternion
    V_l1 = quatrotate(quatinv(l1_quat_total), V_l1_gl');
    omega_l1 = cross(l1_senspoint',V_l1)/norm(l1_senspoint')^2;
    
    %now make quaternion of local delta rotation
    dangle_vec_l1 = omega_l1*dt; 
    dangle_l1 = sqrt(dangle_vec_l1(1)^2 + dangle_vec_l1(2)^2 + dangle_vec_l1(3)^2); 
    
    dqw = cos(dangle_l1/2)/dangle_l1 ;
    dqx = dangle_vec_l1(1) * sin(dangle_l1/2) / dangle_l1;
    dqy = dangle_vec_l1(2) * sin(dangle_l1/2) /dangle_l1;
    dqz = dangle_vec_l1(3) * sin(dangle_l1/2) /dangle_l1;
    
    l1_quat_local = [dqw dqx dqy dqz];
    l1_quat_local;
    l1_quat_total;
    %and now the very beatifule quaterion magic!
    l1_quat_total = quatmultiply (l1_quat_local, l1_quat_total);
    l1_quat_total = quatnormalize(l1_quat_total);
    %according to task - calculate angular velocities (global frame)
    omega_l1_global = cross(l1_senspoint_gl',V_l1_gl)/norm(l1_senspoint_gl')^2;
    omegal1_log = [omegal1_log omega_l1_global];
    
    %and angular accelerations
    epsilon_l1_gl = cross(l1_senspoint_gl', a_l1_gl - cross(omega_l1_global,V_l1_gl))/norm(l1_senspoint_gl')^2;
    epsilon_log = [epsilon_log epsilon_l1_gl];
    
% LINK 2

    [angl2onl1_x, angl2onl1_y, angl2onl1_z]  = quat2angle(l2_quat_onl1_total,'XYZ');
    angl2onl1_xlog = [angl2onl1_xlog angl2onl1_x]; 
    angl2onl1_ylog = [angl2onl1_ylog angl2onl1_y]; 
    angl2onl1_zlog = [angl2onl1_zlog angl2onl1_z]; 
    %to obtain link2 angles in ref to link1 we have to get pure
    %acceleration of link 2 on link 1. For this purpose we have to define
    %start of the link2 (that is link 1 end) acceleration
    %we use formula for acceleration of fixed rotating body - link 1
    
    [l2_base_angx, l2_base_angy, l2_base_angz]  = quat2angle( l2_quat_total_base, 'XYZ');
    l2_base_angles_log = [l2_base_angles_log [l2_base_angx; l2_base_angy; l2_base_angz]];
    R_l2_base = quat2rotm(l2_quat_total_base);
    display(R_l2_base)
    

    
    %we convert data from sensors to the link2-frame
    a1onlink_l2 = a1_l2(step,:)';%if look along Xaxis - this is top sensor
    a2onlink_l2 = rotx(-90)*a2_l2(step,:)';%if look along Xaxis - this is right sensor
    a3onlink_l2 = rotx(-270)*a3_l2(step,:)';%left sensor
    a4onlink_l2 = rotx(-180)*a4_l2(step,:)';%bottom sensor

    a1onlink_l2_gl = quatrotate(l2_quat_total_base, a1onlink_l2');
    a2onlink_l2_gl = quatrotate(l2_quat_total_base, a2onlink_l2');
    a3onlink_l2_gl = quatrotate(l2_quat_total_base, a3onlink_l2');
    a4onlink_l2_gl = quatrotate(l2_quat_total_base, a4onlink_l2');

    a1motion_l2_gl = a1onlink_l2_gl - a_gravity_gl;%this is accelerations with link1 acceleration inside - global acceleration
    a2motion_l2_gl = a2onlink_l2_gl - a_gravity_gl;
    a3motion_l2_gl = a3onlink_l2_gl - a_gravity_gl;
    a4motion_l2_gl = a4onlink_l2_gl - a_gravity_gl;
    
    %and get pseudo-equilibrium acceleration
    a_l2_gl = [0;0;0];    
    a_l2_gl(1) = norm([a1motion_l2_gl(1);a2motion_l2_gl(1);a3motion_l2_gl(1);a4motion_l2_gl(1) ]);
    a_l2_gl(2) = norm([a1motion_l2_gl(2);a2motion_l2_gl(2);a3motion_l2_gl(2);a4motion_l2_gl(2) ]);
    a_l2_gl(3) = norm([a1motion_l2_gl(3);a2motion_l2_gl(3);a3motion_l2_gl(3);a4motion_l2_gl(3) ]);
    
    % to get link 2 motion on link1 we have to substract acceleration of
    % link1, affecting on link 2. For this purpose we have linear
    % acceleration on full rotation formula for sensor point of link2 in
    % global coordinates, so 
    
    %get points of link2 in link1 frame
    l2_senspoint_onl1 = quatrotate(l2_quat_onl1_total, l2_senspoint_zero);
    l2_endpoint_onl1 = quatrotate(l2_quat_onl1_total, l2_endpoint_zero);    
    l2_senspoint_gl = quatrotate(l1_quat_total,l2_senspoint_onl1)+...
                        l1_endpoint_gl;
    l2_endpoint_gl = quatrotate(l1_quat_total,l2_endpoint_onl1)+...
                        l1_endpoint_gl;
                    
                    
    a_l2_l1caused_gl = cross(epsilon_l1_gl, l2_senspoint_gl) + ...
                     +cross( omega_l1_global, ...
                            cross(omega_l1_global, l2_senspoint_gl));
    
                        
    %and pure motion of link2 on(in refernce on) link1 (as if link1 stationary), (gravity
    %acceleration is already considered)

    a_l2_onl1_gl = a_l2_gl - a_l2_l1caused_gl';
    

    %we rotate to link1 reference frame and accumulate velocity of link2 
    %in this frame
    a_l2_refl1 = quatrotate(quatinv(l1_quat_total),a_l2_onl1_gl');
    dV_l2_refl1 = a_l2_refl1 * dt;
    V_l2_refl1 = dV_l2_refl1+ dV_l2_refl1;
    V_l2_onl1 = V_l2_refl1';
    a_l2_onl1 = a_l2_refl1;
    %go to l2 rotated position to obtain angle shift
    V_l2 = quatrotate(quatinv(l2_quat_onl1_total), V_l2_refl1);
    omega_l2 = cross(l2_senspoint_zero',V_l2)/norm(l2_senspoint_zero)^2;


    %integrate angular velocity
    dangle_vec_l2 = omega_l2 * dt;
    dangle_l2 = sqrt(dangle_vec_l2(1)^2 + dangle_vec_l2(2)^2 + dangle_vec_l2(3)^2);
    
    dqw2 = cos(dangle_l2/2)/dangle_l2 ;
    dqx2 = dangle_vec_l2(1) * sin(dangle_l2/2) /dangle_l2;
    dqy2 = dangle_vec_l2(2) * sin(dangle_l2/2) /dangle_l2;
    dqz2 = dangle_vec_l2(3) * sin(dangle_l2/2) /dangle_l2;   
    
    l2_quat_onl1_local = [dqw2 dqx2 dqy2 dqz2];
    l2_quat_onl1_total = quatmultiply(l2_quat_onl1_local, l2_quat_onl1_total);
    l2_quat_onl1_total = quatnormalize(l2_quat_onl1_total);
        
    omega_l2_onl1 = cross(l2_senspoint_onl1',V_l2_onl1)/norm(l2_senspoint_onl1)^2;  
    omega_l2onl1_log = [omega_l2onl1_log omega_l2_onl1];
    
%     epsilon_l2_onl1 = cross(l2_senspoint_onl1',...
%                             a_l2_onl1' - a_l1_end - cross(omega_l2_onl1,V_l2_onl1))...
%                             /norm(l2_senspoint_onl1)^2;

    a_l1_end_gl = cross(epsilon_l1_gl , l1_endpoint_gl') ...
               + cross(omega_l1_global,...
                        cross(omega_l1_global, l1_endpoint_gl'));
    
    a_l1_end_onl1 = quatrotate(quatinv(l1_quat_total), a_l1_end_gl');
                
    epsilon_l2_onl1 = cross(l2_senspoint_onl1',...
                            a_l2_onl1' - a_l1_end_onl1' - cross(omega_l2_onl1,...
                                cross(omega_l2_onl1,l2_senspoint_onl1')))...
                            /norm(l2_senspoint_onl1)^2;
                                              
    epsilon_l2onl1_log = [epsilon_l2onl1_log epsilon_l2_onl1];

        
    %and total rotation of link2 in base coordinates is required for
    %motion acceleration extracting   
    l2_quat_total_base = quatmultiply(l2_quat_onl1_total, l1_quat_total);

    h2 = figure(2);
    xlim([-1.6 1.6])
    ylim([-1.6 1.6])
    zlim([-1.6 1.6])
    grid on
    axis equal
    %black rod - link1 body
    plot3([0; l1_endpoint_gl(1); l2_endpoint_gl(1)],...
        [0; l1_endpoint_gl(2); l2_endpoint_gl(2)],...
        [0; l1_endpoint_gl(3); l2_endpoint_gl(3)],...
        'color','black','linewidth',5)
    hold on
    %arrows of accelerations
%     mArrow3(l1_senspoint_gl',a1motion_l1_gl/10,'color','red'); hold on
%     mArrow3(l1_senspoint_gl',a2motion_l1_gl/10,'color','blue'); hold on
%     mArrow3(l1_senspoint_gl',a3motion_l1_gl/10,'color','magenta'); hold on
%     mArrow3(l1_senspoint_gl',a4motion_l1_gl/10,'color','cyan'); hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    % redball - base
    scatter([0],[0], 55,'red','filled')
    hold on
    scatter3([l1_endpoint_gl(1)],[l1_endpoint_gl(2)],[l1_endpoint_gl(3)], 55,'green','filled')

%     title('Link1 initial pose with MOTION accelerations /10.0(gravity excluded)')
    title('Robot poses')

    grid on
    axis equal
    xlim([-0.6 0.6])
    ylim([-0.6 0.6])
    zlim([-0.6 0.6])
%     hold off

%     sqrt(...
%     (l1_endpoint_gl(1) - l2_endpoint_gl(1))^2+...
%     (l1_endpoint_gl(2) - l2_endpoint_gl(2))^2+...
%     (l1_endpoint_gl(3) - l2_endpoint_gl(3))^2)


    
end


%plot for Link1
%-------------------------------------
figure(10)
% title('SRD1 file. Euler angles (XYZ sequence) for link 1 in reference to base.')

subplot(3,1,1)
plot(rad2deg(angl1_xlog), 'bo')
xlabel('Measurment Number')
ylabel('X angle, degrees')
grid on 

subplot(3,1,2)
plot(rad2deg(angl1_ylog), 'bo')
xlabel('Measurment Number')
ylabel('Y angle, degrees')
grid on 

subplot(3,1,3)
plot(rad2deg(angl1_zlog),'bo')
xlabel('Measurment Number')
ylabel('Z angle, degrees')
grid on 

set(gcf,'NextPlot','add');
axes;
h = title('Euler angles (XYZ sequence) for link 1 in reference to base.');
set(gca,'Visible','off');
set(h,'Visible','on');


figure(11)

subplot(3,1,1)
%plot(rad2deg(omega_log(1,:)), 'bo')
plot((omegal1_log(1,:)/(2*pi)), 'bo')
xlabel('Measurment Number')
%ylabel('X angular velocity, degrees/sec')
ylabel('X angular velocity, revolution/sec')
grid on 

subplot(3,1,2)
% plot(rad2deg(omega_log(2,:)), 'bo')
plot((omegal1_log(2,:)/(2*pi)), 'bo')
xlabel('Measurment Number')
% ylabel('Y angular velocity, degrees/sec')
ylabel('Y angular velocity, revolution/sec')
grid on 

subplot(3,1,3)
% plot(rad2deg(omega_log(3,:)),'bo')
plot(omegal1_log(3,:)/(2*pi),'bo')
xlabel('Measurment Number')
% ylabel('Z angular velocity, degrees/sec')
ylabel('Z angular velocity, revolution/sec')
grid on 

set(gcf,'NextPlot','add');
axes;
h2 = title('Angular velocitiy for link 1');
set(gca,'Visible','off');
set(h2,'Visible','on');


figure(12)

subplot(3,1,1)
% plot(rad2deg(epsilon_log(1,:)), 'bo')
plot(epsilon_log(1,:)/(2*pi), 'bo')
xlabel('Measurment Number')
ylabel('X angular acceleration, revolution/sec^2')
grid on 

subplot(3,1,2)
%plot(rad2deg(epsilon_log(2,:)), 'bo')
plot(epsilon_log(2,:)/(2*pi), 'bo')
xlabel('Measurment Number')
%ylabel('Y angular acceleration, degrees/sec^2')
ylabel('Y angular acceleration, revolutions/sec^2')
grid on 

subplot(3,1,3)
% plot(rad2deg(epsilon_log(3,:)),'bo')
plot(epsilon_log(3,:)/(2*pi),'bo')
xlabel('Measurment Number')
% ylabel('Z angular acceleration, degrees/sec^2')
ylabel('Z angular acceleration, revolutions/sec^2')
grid on 

set(gcf,'NextPlot','add');
axes;
h3 = title('Angular acceleration for link 1');
set(gca,'Visible','off');
set(h3,'Visible','on');

%--------------------------------------------------------------
%Plots for link 2

figure(20)

subplot(3,1,1)
plot(rad2deg(angl2onl1_xlog), 'bo')
xlabel('Measurment Number')
ylabel('X angle, degrees')
grid on 

subplot(3,1,2)
plot(rad2deg(angl2onl1_ylog), 'bo')
xlabel('Measurment Number')
ylabel('Y angle, degrees')
grid on 

subplot(3,1,3)
plot(rad2deg(angl2onl1_zlog),'bo')
xlabel('Measurment Number')
ylabel('Z angle, degrees')
grid on 

set(gcf,'NextPlot','add');
axes;
h = title('Euler angles (XYZ sequence) for link 2 in reference to link1.');
set(gca,'Visible','off');
set(h,'Visible','on');


figure(21)

subplot(3,1,1)
%plot(rad2deg(omega_log(1,:)), 'bo')
plot((omega_l2onl1_log(1,:)/(2*pi)), 'bo')
xlabel('Measurment Number')
%ylabel('X angular velocity, degrees/sec')
ylabel('X angular velocity, revolution/sec')
grid on 

subplot(3,1,2)
% plot(rad2deg(omega_log(2,:)), 'bo')
plot((omega_l2onl1_log(2,:)/(2*pi)), 'bo')
xlabel('Measurment Number')
% ylabel('Y angular velocity, degrees/sec')
ylabel('Y angular velocity, revolution/sec')
grid on 

subplot(3,1,3)
% plot(rad2deg(omega_log(3,:)),'bo')
plot(omega_l2onl1_log(3,:)/(2*pi),'bo')
xlabel('Measurment Number')
% ylabel('Z angular velocity, degrees/sec')
ylabel('Z angular velocity, revolution/sec')
grid on 

set(gcf,'NextPlot','add');
axes;
h2 = title('Angular velocitiy for link 2 in reference to link1');
set(gca,'Visible','off');
set(h2,'Visible','on');


figure(22)

subplot(3,1,1)
% plot(rad2deg(epsilon_log(1,:)), 'bo')
plot(epsilon_l2onl1_log(1,:)/(2*pi), 'bo')
xlabel('Measurment Number')
ylabel('X angular acceleration, revolution/sec^2')
grid on 

subplot(3,1,2)
%plot(rad2deg(epsilon_log(2,:)), 'bo')
plot(epsilon_l2onl1_log(2,:)/(2*pi), 'bo')
xlabel('Measurment Number')
%ylabel('Y angular acceleration, degrees/sec^2')
ylabel('Y angular acceleration, revolutions/sec^2')
grid on 

subplot(3,1,3)
% plot(rad2deg(epsilon_log(3,:)),'bo')
plot(epsilon_l2onl1_log(3,:)/(2*pi),'bo')
xlabel('Measurment Number')
% ylabel('Z angular acceleration, degrees/sec^2')
ylabel('Z angular acceleration, revolutions/sec^2')
grid on 

set(gcf,'NextPlot','add');
axes;
h3 = title('Angular acceleration of link 2 in reference to link 1');
set(gca,'Visible','off');
set(h3,'Visible','on');


%------------------LINK 2 BASE FRAME-------------
figure(23)

subplot(3,1,1)
plot(rad2deg(l2_base_angles_log(1,:)), 'bo')
xlabel('Measurment Number')
ylabel('X angle, degrees')
grid on 

subplot(3,1,2)
plot(rad2deg(l2_base_angles_log(2,:)), 'bo')
xlabel('Measurment Number')
ylabel('Y angle, degrees')
grid on 

subplot(3,1,3)
plot(rad2deg(l2_base_angles_log(3,:)),'bo')
xlabel('Measurment Number')
ylabel('Z angle, degrees')
grid on 

set(gcf,'NextPlot','add');
axes;
h = title('Euler angles (XYZ sequence) for link 2 in reference to base');
set(gca,'Visible','off');
set(h,'Visible','on');

