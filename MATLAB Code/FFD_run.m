%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for testing FFD and methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc

% instantiate a new FFD object
FFD_test = FFD();

% testing iterateUrStar method in FFD
UstarTest1 = FFD_test.Ustar;
FFD_test.computeUStar;
UstarTest2 = FFD_test.Ustar;
ArStarTest = full(FFD_test.ArStar);
NrTest = reshape(FFD_test.Nr, length(FFD_test.rbar), ...
          length(FFD_test.zbar))';
AzStarTest = full(FFD_test.AzStar);
NzTest = reshape(FFD_test.Nz, length(FFD_test.rbar), ...
          length(FFD_test.zbar))';
