function Task_4_5()
    clc
    clear all
    
    s_le = [-0.1878 -0.1021];
    s_te = [-0.0262 0.0147];
    a_le = [0.0002 -0.0024];
    a_te = [1 0];
    f_le = [1.0397 -0.0207];
    f_te = [1.2666 -0.1257];
    
    r_lr = 0.022;
    omega = 0.2;
    
    t_start = tic;
    tic
    ref_local('Mesh_1_26_2023.gri',s_le(1),s_le(2),r_lr+0.015,omega,'Mesh_1_26_2023_lr_sle.gri')
    ref_local('Mesh_1_26_2023_lr_sle.gri',s_te(1),s_te(2),r_lr,omega,'Mesh_1_26_2023_lr_sa.gri')
    ref_local('Mesh_1_26_2023_lr_sa.gri',a_le(1),a_le(2),r_lr,omega,'Mesh_1_26_2023_lr_sa_ale.gri')
    ref_local('Mesh_1_26_2023_lr_sa_ale.gri',a_te(1),a_te(2),r_lr,omega,'Mesh_1_26_2023_lr_sa_aa.gri')
    ref_local('Mesh_1_26_2023_lr_sa_aa.gri',f_le(1),f_le(2),r_lr,omega,'Mesh_1_26_2023_lr_sa_aa_fle.gri')
    ref_local('Mesh_1_26_2023_lr_sa_aa_fle.gri',f_te(1),f_te(2),r_lr+0.015,omega,'Mesh_1_26_2023_lr.gri')
    
    ref_uniform('Mesh_1_26_2023_lr.gri','Mesh_1_26_2023_rf1.gri')
    ref_uniform('Mesh_1_26_2023_rf1.gri','Mesh_1_26_2023_rf2.gri')
    ref_uniform('Mesh_1_26_2023_rf2.gri','Mesh_1_26_2023_rf3.gri')
    t_stop = toc;
    
    disp ('ho gaya');
end