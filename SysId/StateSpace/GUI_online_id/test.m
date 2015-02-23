clc;
clear;
close all;

load('debugValues');

weights = [1,3,1,35];

sim_MpcVsLqr(realModel, lqrModel, mpcModel, predHor_steps, ...
                      weights, deltas, constraints, typeTack);