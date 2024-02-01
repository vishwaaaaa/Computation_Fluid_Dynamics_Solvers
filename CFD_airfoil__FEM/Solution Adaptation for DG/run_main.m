clear all
clc

%solnAdaptation_f_DG(Input .gri file,solution order,Input state file,Error Metric at edges file,Output .gri file,Output state file)

solnAdaptation_f_DG('c0.gri',0,'c0_p0_q1_U.txt','adaptMetric_q1_p0.txt','c0_p0_adapt.gri','c0_p0_state_adapt.txt')
solnAdaptation_f_DG('c0.gri',1,'c0_p1_q1_U.txt','adaptMetric_q1_p1.txt','c0_p1_adapt.gri','c0_p1_state_adapt.txt')