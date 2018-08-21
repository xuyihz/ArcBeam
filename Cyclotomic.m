%% Cyclotomic
% 
% 
% Xu Yi, 2018

%%
close all; clear; clc;

%%
lamda = 0.995;      % 以直代曲 直线段长度和与圆弧长度之比
                    % 直线段长度 sin(theta_i/2)*R*2
                    % 圆弧长度   theta*R
theta = pi/2;       % 1/4圆
n = 1;              % 初始化 角度份数
theta_i = theta;    % 初始化 单位角度
while sin(theta_i/2)*2/theta_i < lamda
    n = n+1;
    theta_i = theta/n;
end
disp(n);