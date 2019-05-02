%Initial position and angular co-ordinates of the robot 
x0 = 1; 
y0 = 1;
theta0 = 0*pi/180;

%Final position and angular co-ordinates of the robot 
xf = 20;
yf = 20;
thetaf = 45*pi/180;

%Initial and final velocity co-ordinates 
xdot0 = 0;
xdotf = 5;
k_t0_dot = 0;
k_tf_dot = 10;

%Initial and final time 
tf = 10;
t0 = 2;

touf= 3;

%Order of the Brienstien Polynomial
n = 3;

%Brienstien Co-eff for t0 and tf
B0_t0 = 1;
B1_t0 = 0;
B2_t0 = 0;
B3_t0 = 0;

B0_tf = 0;
B1_tf = 0;
B2_tf = 0;
B3_tf = 1;

B0_t0_dot = -3/8;
B1_t0_dot = 3/8;
B2_t0_dot = 0;
B3_t0_dot = 0;

B0_tf_dot = 0;
B1_tf_dot = 0;
B2_tf_dot = -3/8;
B3_tf_dot = 3/8;



W_x0 = x0;
W_x3 = xf;

A = [xdot0-W_x0*B0_t0_dot- W_x3*B3_t0_dot ; xdotf-W_x0*B0_tf_dot- W_x3*B3_tf_dot];
B = [B1_t0_dot  B2_t0_dot; B1_tf_dot  B2_tf_dot ];

W = pinv(B)*A; %weights matrix for x(t)

k0 = tan(theta0);
kf = tan(thetaf);

W_k0 = k0;
W_k3 = kf;

C = [k_t0_dot-W_k0*B0_t0_dot-W_k3*B3_t0_dot ; k_tf_dot-W_k0*B0_tf_dot-W_k3*B3_tf_dot ];

W_k = pinv(B)*C; %weights matrix for k(t)
 
syms t1

mu1 = (t1-t0)/(tf-t0);

B0_dot_mu = nchoosek(3,0)*((-3*((1-mu1)^2))/(tf-t0));
B1_dot_mu = nchoosek(3,1)*(((-2*mu1*(1-mu1))+((1-mu1)^2))/(tf-t0));
B2_dot_mu = nchoosek(3,2)*(((-mu1^2)+(2*mu1*(1-mu1)))/(tf-t0));
B3_dot_mu = nchoosek(3,3)*((3*mu1^2)/(tf-t0));

B0_mu = nchoosek(3,0)*((1-mu1)^3);
B1_mu = nchoosek(3,1)*((1-mu1)^2)*mu1;
B2_mu = nchoosek(3,2)*((1-mu1)^1)*(mu1^2);
B3_mu = nchoosek(3,3)*((1-mu1)^0)*(mu1^3);


F_0 = (B0_dot_mu*W_x0*B0_mu)+(B1_dot_mu*W(1,1)*B0_mu)+(B2_dot_mu*W(2,1)*B0_mu)+(B3_dot_mu*W_x3*B0_mu);
F_1 = (B0_dot_mu*W_x0*B1_mu)+(B1_dot_mu*W(1,1)*B1_mu)+(B2_dot_mu*W(2,1)*B1_mu)+(B3_dot_mu*W_x3*B1_mu);
F_2 = (B0_dot_mu*W_x0*B2_mu)+(B1_dot_mu*W(1,1)*B2_mu)+(B2_dot_mu*W(2,1)*B2_mu)+(B3_dot_mu*W_x3*B2_mu);
F_3 = (B0_dot_mu*W_x0*B3_mu)+(B1_dot_mu*W(1,1)*B3_mu)+(B2_dot_mu*W(2,1)*B3_mu)+(B3_dot_mu*W_x3*B3_mu);

F_0_t0 = double(subs(int(F_0),t0));
F_1_t0 = double(subs(int(F_1),t0));
F_2_t0 = double(subs(int(F_2),t0));
F_3_t0 = double(subs(int(F_3),t0));

F_0_tf = double(subs(int(F_0),tf));
F_1_tf = double(subs(int(F_1),tf));
F_2_tf = double(subs(int(F_2),tf));
F_3_tf = double(subs(int(F_3),tf));





D = [y0-(W_k0*F_0_t0)-(W_k3*F_3_t0); yf-(W_k0*F_0_tf)-(W_k3*F_3_tf)];
E = [F_1_t0  F_2_t0 ; F_1_tf  F_2_tf];
W_y = pinv(E)*D; %weights for y(t)

x_t=[];
k_t=[];
y_t=[];
theta_t =[];
B_r = [];
x_t_vel = [];
y_t_vel = [];
x_t_vel_new = [];
y_t_vel_new = [];
omega =[];
a = [];


   F0_int = int(F_0);
   F1_int = int(F_1);
   F2_int = int(F_2);
   F3_int = int(F_3);

 x_t_eq = W_x0*B0_mu + W(1,1)*B1_mu + W(2,1)*B2_mu + W_x3*B3_mu;
 y_t_eq = y0+W_k0*F0_int + W_y(1,1)*F1_int + W_y(2,1)*F2_int + W_k3*F3_int;
 k_t_eq =W_k0*B0_mu + W_k(1,1)*B1_mu+ W_k(2,1)*B2_mu + W_k3*B3_mu;
 theta_t_eq = atan(k_t_eq);
 
 x_t_velocity = diff(x_t_eq);
 y_t_velocity = diff(y_t_eq);
 theta_t_velocity = diff(theta_t_eq);

 %time loop between t0 and tf for running and plotting the brienstien
 %polynomial based trajectory planning and generation


for t = 2:0.01:10
    
    for i = 0:n
        mu =  (t-t0)/(tf-t0);
        f = (1-mu)^(n-i);
        g = (mu)^(i);
        B_r= [B_r nchoosek(n,i)*f*g];
       
    end
    
    
  F0 = double(int(F_0,t0,t));
  F1 = double(int(F_1,t0,t));
  F2 = double(int(F_2,t0,t));
  F3 = double(int(F_3,t0,t));
  
   
   
  % x(t) y(t) theta(t)   
  x_t = [x_t W_x0*B_r(end-3) + W(1,1)*B_r(end-2) + W(2,1)*B_r(end-1) + W_x3*B_r(end)];
  k_t = [k_t  W_k0*B_r(end-3) + W_k(1,1)*B_r(end-2)+ W_k(2,1)*B_r(end-1) + W_k3*B_r(end)]; 
  y_t = [y_t y0+W_k0*F0 + W_y(1,1)*F1 + W_y(2,1)*F2 + W_k3*F3 ];
  theta_t =  [theta_t atan(k_t(end))];
  x_t_vel = [x_t_vel double(subs(x_t_velocity,t))];
  y_t_vel = [y_t_vel double(subs(y_t_velocity,t))];
  
  x_t_vel_i = double(subs(x_t_velocity,t));
  y_t_vel_i =  double(subs(y_t_velocity,t));
  
  tc = 2;
  
  a = (tf-tc)/(touf-tc);
  x_t_vel_new = [x_t_vel_new a*x_t_vel_i];
  y_t_vel_new = [y_t_vel_new a*y_t_vel_i];
  omega = [omega double(subs(theta_t_velocity,t))];
  
end


% %Plotting of the trajectory 
% plot(x_t,y_t)
% hold on
% p = plot(x_t(1),y_t(1),'o','MarkerFaceColor','red','MarkerSize',8);
% hold off
% 
% for k = 2:length(x_t)
%     p.XData = x_t(k);
%     p.YData = y_t(k);
%     drawnow 
% end

%% plotting
x_t_new = [];
y_t_new = [];
theta0 = theta_t(1);
for l = 2:length(x_t)
    d =  (x_t_vel_new(l-1)*(0.01));
    x_t_new = [x_t_new x_t(l-1) + d];
    y_t_new = [y_t_new  y_t(l-1)+(y_t_vel_new(l-1)*(0.01))];
  
    axis([0 25 0 25]); 
    hold on 
    plot(x_t_new,y_t_new,'b');
    plot(x_t_new(l-1),y_t_new(l-1),'o','MarkerFaceColor','red','MarkerSize',8);
    plot(x_t(l-1),y_t(l-1),'o','MarkerFaceColor','g','MarkerSize',8);
    legend('Bernstien path','Time scaled', 'Non-Time Scaled');
    pause(0.01);
    cla
end
hold off
plot(x_t_vel_new);
hold on
plot(x_t_vel);
legend('Time scaled Velocity', 'Non-Time Scaled Velocity','Location','northwest');
