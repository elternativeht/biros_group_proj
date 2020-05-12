%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for testing FFD and methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc

% instantiate a new FFD object using
FFD_test = FFD();
% this will initialize FFD_test with all properties and methods defined in
% FFD with all default definitions. it will also run the object 
% initialization method.

% run any of the methods externaly using
var1 = 1; var2 = 2; 
x1 = FFD_test.exampleFunction(var1, var2);
% or 
x2 = exampleFunction(FFD_test, var1, var2);
% for functions that have outputs and
rbar_old = FFD_test.rbar; zbar_old = FFD_test.zbar;
FFD_test.exampleOperation(var1, var2);
rbar_new1 = FFD_test.rbar; zbar_new1 = FFD_test.zbar;
% or
exampleOperation(FFD_test, var1, var2);
rbar_new2 = FFD_test.rbar; zbar_new2 = FFD_test.zbar;
% for operations. in an operation, all associated properties will be
% updated. e.g. rbar and zbar in this case

% testing stress computation
n = length(FFD_test.zbar);
m = length(FFD_test.rbar);
nm = n*m;
ur = logspace(0, 1, nm)';
uz = ones(nm, 1);
[Srr, Szz, tauR, tauZ] = FFD_test.stress(ur, uz);

% testing patchVelocity
% FFD_test.Urbar = 0.5*ones(FFD_test.zMaxIndex+1, FFD_test.rMaxIndex);
% FFD_test.patchVelocity(2);
uTest = full(FFD_test.u);
vTest = full(FFD_test.v);

% testing captureSpeed and captureStress
FFD_test.captureSpeed(2);
FFD_test.captureStress(2);

% testing animateSpeed and animateStress
FFD_test.animateSpeed();
FFD_test.animateStress();

% testing animateuMag
% FFD_test.animateUMag();

% testing computeNr and computeNz
% computeNr(FFD_test);
% computeNz(FFD_test);
% 
% test1 = reshape(FFD_test.Nr, size(FFD_test.Urbar, 2), ...
%                          size(FFD_test.Urbar, 1))';
% test2 = reshape(FFD_test.Nz, size(FFD_test.Uzbar, 2), ...
%                size(FFD_test.Uzbar, 1))';
           
% testing computeDstar
% computeDstar(FFD_test);
% test3 = reshape(FFD_test.Dstar, size(FFD_test.Pbar, 2), ...
%                          size(FFD_test.Pbar, 1))';

% testing iterateUrStar method in FFD
% UstarTest1 = FFD_test.Ustar;
% FFD_test.computeUStar;
% UstarTest2 = FFD_test.Ustar;
% ArStarTest = full(FFD_test.ArStar);
% NrTest = reshape(FFD_test.Nr, length(FFD_test.rbar), ...
%           length(FFD_test.zbar))';
% AzStarTest = full(FFD_test.AzStar);
% NzTest = reshape(FFD_test.Nz, length(FFD_test.rbar), ...
%           length(FFD_test.zbar))';
