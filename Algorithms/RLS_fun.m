function [w,Sd,updates] = RLS_fun(w,Sd,z,lambda,e)
%Signals, Multimedia and Telecommunications Lab
%
%This function updates the weights using the RLS algorithm
%
%Author: Felipe Barboza da Silva       felipe.silva@smt.ufrj.br

updates = 1;
psi = Sd*z;
Sd = (1/lambda)*(Sd-((psi*psi')/(lambda+psi'*z)));
w = w + e*Sd*z;