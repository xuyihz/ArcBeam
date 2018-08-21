%% Cyclotomic
% 
% 
% Xu Yi, 2018

%%
close all; clear; clc;

%%
lamda = 0.995;      % ��ֱ���� ֱ�߶γ��Ⱥ���Բ������֮��
                    % ֱ�߶γ��� sin(theta_i/2)*R*2
                    % Բ������   theta*R
theta = pi/2;       % 1/4Բ
n = 1;              % ��ʼ�� �Ƕȷ���
theta_i = theta;    % ��ʼ�� ��λ�Ƕ�
while sin(theta_i/2)*2/theta_i < lamda
    n = n+1;
    theta_i = theta/n;
end
disp(n);