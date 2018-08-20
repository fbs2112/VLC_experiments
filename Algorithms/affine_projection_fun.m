function [w,updates] = affine_projection_fun(w,z,muAP,e)
%Signals, Multimedia and Telecommunications Lab
%
%This function updates the weights using the Affine Projection algorithm 
%
%Author: Felipe Barboza da Silva       felipe.silva@smt.ufrj.br

updates = 1;
gamma = 1e-8;
reuseWindowLength = size(z,2) - 1;

 w = w + muAP*z*((z'*z+gamma*eye(reuseWindowLength+1))\eye(reuseWindowLength+1))*e;