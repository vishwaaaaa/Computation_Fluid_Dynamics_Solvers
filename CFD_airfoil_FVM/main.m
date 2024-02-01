clc
clear all

gri_in_file ='c0_ED_02_It4.gri';                         % input gri file
%fstatein_1 = "C:\Users\benba\Documents\C++\AERO623Proj2\FirstOrderAdapt4SolRaw.txt";  % 1st order solution on input file
fstatein_2 = "C:\Users\benba\Documents\C++\AERO623Proj2\SecondOrderAdapt4SolRaw.txt"; % 2nd order solution in input file
gamma = 1.4;                                      % fluid property; no need to change it
%gri_out_file_1 = 'c0_ED_01_It5.gri';              % gri file refined based on 1st order solution
gri_out_file_2 = 'c0_ED_02_It5.gri';              % gri file refined based on 2nd order solution
%fstateout_1 = 'C:\Users\benba\Documents\C++\AERO623Proj2\FirstOrderAdapt4Sol.txt';  % 1st order solution mapped on gri file refined based on 1st order solution
fstateout_2 = 'C:\Users\benba\Documents\C++\AERO623Proj2\SecondOrderAdapt4Sol.txt'; % 2nd order solution mapped on gri file refined based on 2nd order solution

%solnAdaptation(gri_in_file,fstatein_1,gamma,gri_out_file_1,fstateout_1);
%disp ('Mesh Adaptation-First Order-Done');
solnAdaptation(gri_in_file,fstatein_2,gamma,gri_out_file_2,fstateout_2);
disp ('Mesh Adaptation-Second Order-Done');