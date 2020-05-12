%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for testing FFD and methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc
delete all
% % instantiate a new FFD object using
% FFD_test = FFD();
% % this will initialize FFD_test with all properties and methods defined in
% % FFD with all default definitions. it will also run the object 
% % initialization method.
% 
% % run any of the methods externaly using
% var1 = 1; var2 = 2; 
% x1 = FFD_test.exampleFunction(var1, var2);
% % or 
% x2 = exampleFunction(FFD_test, var1, var2);
% % for functions that have outputs and
% rbar_old = FFD_test.rbar; zbar_old = FFD_test.zbar;
% FFD_test.exampleOperation(var1, var2);
% rbar_new1 = FFD_test.rbar; zbar_new1 = FFD_test.zbar;
% % or
% exampleOperation(FFD_test, var1, var2);
% rbar_new2 = FFD_test.rbar; zbar_new2 = FFD_test.zbar;
% % for operations. in an operation, all associated properties will be
% % updated. e.g. rbar and zbar in this case
% 
% % testing iterateUrStar method in FFD
% 
% UstarTest1 = FFD_test.Ustar;
% FFD_test.computeUStar;
% UstarTest2 = FFD_test.Ustar;
% ArStarTest = full(FFD_test.ArStar);
% NrTest = reshape(FFD_test.Nr, length(FFD_test.rbar), ...
%           length(FFD_test.zbar))';
% AzStarTest = full(FFD_test.AzStar);
% NzTest = reshape(FFD_test.Nz, length(FFD_test.rbar), ...
%           length(FFD_test.zbar))';


FFD_test = FFD(); FFD_test.computeUStar();
FFD_test.computeu();


      
      
      
      
      
      