function [minus] = ATT_EKF_DOT(pre)

% Yazan Chihabi, 2018

q_pre = pre(1:4);
bias_pre = pre(5:7);
w_pre = pre(8:10);


wx = w_pre(1); wy = w_pre(2); wz = w_pre(3);

w_skew = [0 -wz wy; wz 0 -wx; -wy wx 0];

Omega = [-w_skew w_pre; -w_pre' 0];


q_dot = 0.5*Omega*q_pre;



tau = 300;
bias_dot = -(1/tau)*bias_pre; 
minus = [q_dot; bias_dot];

end