%% Sopravvivenza_AIRO! 
% Author: Massimo, Leonardo, Paolo, Francesco
% This script take as input the direct kinemantics (task kinematics) 
% of a robot manipulator and compute the regressor matrix and 
% the new DH parameters

clc
clear all
close all

syms da dalpha dd dtheta theta a alpha d l1 l2 q1 q2 alpha1 alpha2 d1 d2 zero 

%% DIRECT KINEMATICS
% Write you task kinematics

r = [l1*cos(q1) + l2*cos(q1+q2); 
     l1*sin(q1) + l2*sin(q1+q2)] 

%% COMPUTATION OF REGRESSOR MATRIX

% Put "zero" if you don't want calibrate that DH parameters
% Otherwise, write symbolically the DH parameters 
% (EX. "l1, q1..." AND SO ON!)

alpha = [zero;
         zero];

a = [l1;
     l2];

d = [zero;
     zero];

theta = [zero;
         zero];

% REMEMBER: comment the appropriate derivatives based on your DH
% calibration parameters

%% SINGLE EXPERIMENT

jac_alpha = jacobian(r,alpha);
jac_a = jacobian(r,a);
jac_d = jacobian(r,d);
jac_theta = jacobian(r,theta);

PHI_sym = [];
PHI_full_val = [];
PHI_full_sym = [];

if not(jac_alpha == zeros(2,2)) 
    PHI_sym = [PHI_sym,jac_alpha];
end

if not(jac_a == zeros(2,2))
    PHI_sym = [PHI_sym,jac_a];
end

if not(jac_d == zeros(2,2))
    PHI_sym = [PHI_sym,jac_d];
end

if not(jac_theta == zeros(2,2))
    PHI_sym = [PHI_sym, jac_theta];
end

% We have wrote the symbolically form of REGRESSOR MATRIX (1 EXPERIMENT!)
fprintf("This is the regressor matrix: \n")

PHI_sym

%% MULTIPLE EXPERIMENT

% COMPUTATION OF REQUIRED MATRICIES (FULL EVALUATE MEASURES AND REGRESSOR)
% These are symbolic parameters used for evaluation of FULL MATRICIES
% ATTENTION: CHANGE THESE IN BASE OF EVALUATION TASK 
% !!! RESPECT ORDER -> alpha,a,d,theta !!! 
% EX. You want to evaluate a, d? Ok, write -> phi_sym=[a,d]
phi_sym = [l1,l2,q1,q2];

% Change p_meas, q_meas and phi_nom depending on your exercise
p_meas = [2;0;0;2;1.6925;0.7425;1.7218;0.6718];
% !!! ATTENTION ORDER l1,l2,q1,q2 !!!
phi_val = {[1,1,0,0], [1,1,pi/2,0], [1,1,pi/4,-pi/4], [1,1,0,pi/4]};
phi_nom = {[1; 1]};  % phi_nom = {[alpha],[a],[d],[theta]}
Dr_full = [];
l = 4; % Number of measurements

PHI_full_val = subs_airo(PHI_sym, phi_sym, phi_val, l);
r_sub = subs_airo(r, phi_sym, phi_val, l);
Dr_full = p_meas - r_sub;

PHI_full_val = round(vpa(PHI_full_val),3)
Dr_full = round(vpa(Dr_full),6)

% ITERATION FOR DELTA PHI

eps = 10^-6;
Dphi = eps + 1;
% THESE ARE THE COSTANT (NOT CHANGE AFTER CALIBRATION) VALUES
phi_val = {[0,0], [pi/2,0], [pi/4,-pi/4], [0,pi/4]};
i=0;
condition=1;
old_max = 10^10
phi_prime=phi_nom;
while (abs(max(Dphi)) > eps)
    i=i+1

    pinv_PHI = round(vpa(pinv(PHI_full_val)),3)
    Dphi = round(vpa(pinv_PHI*Dr_full),6)
    fprintf("These are the calibrated DH parameters at first step: \n")
    phi_prime = phi_prime+Dphi
    PHI_sym_tmp = subs(PHI_sym, a, phi_prime) % Evaluate l1,l2
    PHI_full_val = subs_airo(PHI_sym_tmp, [q1,q2], phi_val, l) % Evaluate theta
    
    r_sub_tmp = subs(r, a, phi_prime); % Evaluate l1,l2
    r_sub = subs_airo(r_sub_tmp, [q1,q2], phi_val, l); % Evaluate q1,q2
    Dr_full = p_meas - r_sub;
    Dr_full = round(vpa(Dr_full),6)
    
end


%% FUNCTION TO DO SUBSTITUTION

function subs_mat = subs_airo (mat_sym, symb, val, index)
    mat_out = [];
    for i=1:index
        new_row = subs(mat_sym, symb, val{i});
        mat_out = [mat_out; new_row];
    end
    subs_mat = mat_out;
end





