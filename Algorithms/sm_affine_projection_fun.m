function [w,updates] = sm_affine_projection_fun(w,z,e,threshold)
%Signals, Multimedia and Telecommunications Lab
%
%This function updates the weights using the set-membership Affine 
%Projection algorithm
%
%Author: Felipe Barboza da Silva       felipe.silva@smt.ufrj.br

gamma = 1e-8;
reuseWindowLength = size(z,2) - 1;
u = zeros(reuseWindowLength + 1,1);
u(1) = 1;

error = abs(e(1));
updates = 0;
if error > threshold
    muAP = 1 - threshold/error; 
    w = w + muAP*z*((z'*z+gamma*eye(reuseWindowLength+1))\eye(reuseWindowLength+1))*u*e;
    updates = 1;
end