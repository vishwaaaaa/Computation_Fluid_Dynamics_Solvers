%EcoLine
%AEROSP 623 
%Project 1
clear all; close all; clc;
%Main Driver File

%Task 1
A = load('GmshOutput.mat');
makegri(A.N,A.E);
clear A

%Tasks 2 and 3
[I2E, B2E, In, Bn, A, Ee] = Task2_Task3('Mesh_1_26_2023.gri');

%Tasks 4 and 5
Task_4_5();
[V, E2N,~] = E2N_func('Mesh_1_26_2023_rf1.gri');
[Ee1,~] = verify(E2N, V);
[V, E2N,~] = E2N_func('Mesh_1_26_2023_rf2.gri');
[Ee2,~] = verify(E2N, V);
[V, E2N,~] = E2N_func('Mesh_1_26_2023_rf3.gri');
[Ee3,~] = verify(E2N, V);
