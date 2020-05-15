%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kaden Plewe
% Jinghu Hu
% Aubrey McCutchan
% Spring 2020 ME 382 N Semester Project
% Prof. George Biros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates an object to represent a particle storage bin system
% that will be hydraulically simulated to produce a velocity and stress
% field that represents a particle discharge process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

% create funnel flow discharge object
binFFD = FFD();

% set desired geometric/fluid parameters for the system
binFFD.H = 1;               % (m) height of the storage bin
binFFD.a0 = 0.04;           % non-dimensional diameter of bin bottom opening
binFFD.b = 0.4;             % non-dimensional outer diameter of bin
binFFD.Uinf = 1e-5;         % (m/s) free-stream velocity at top of bin

% set mesh lengths
binFFD.drbar = 0.001;       % non-dimensional radial mesh length
binFFD.dzbar = 0.05;        % non-dimensional z mesh length

% use von Neumann analysis to select sufficient time step
binFFD.dtau = 1e-7;         % non-dimensional simulation time step
binFFD.tauEnd = 1e-6;       % non-dimensional simulation end time

% reinitialize binFFD to account for variable changes
reInitObj(binFFD);

% run funnel flow discharge simulation
simulateFFD(binFFD);










