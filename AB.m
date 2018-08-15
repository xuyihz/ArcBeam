%% Generate MGT file
% main M file
% 
% Xu Yi, 2018

%%
close all; clear; clc;

%% 
fileID = fopen('ArcBeam.mgt','w');   % Open or create new file for writing. Discard existing contents, if any.
addpath(genpath('coor_fun'))    % ����·���м���coor_fun�ļ��м����������ļ���

%% append initial conditions
AB_init(fileID);

%% append model file
AB_model(fileID);

%%
fprintf(fileID,'*ENDDATA');

%%
fclose('all');
