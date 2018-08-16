%% Generate MGT file
% main M file
% 
% Xu Yi, 2018

%%
close all; clear; clc;

%% 
fileID = fopen('CurvedBeam.mgt','w');   % Open or create new file for writing. Discard existing contents, if any.
addpath(genpath('coor_fun'))    % ����·���м���coor_fun�ļ��м����������ļ���

%% append initial conditions
CB_init(fileID);

%% append model file
CB_model(fileID);

%%
fprintf(fileID,'*ENDDATA');

%%
fclose('all');
