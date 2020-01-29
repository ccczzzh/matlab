%% Return Alpha Channel of PNG Image  
% Return the alpha channel of the sample image, |peppers.png|.   

% Copyright 2015 The MathWorks, Inc.


%%  
[X,map,alpha] = imread('peppers.png');
alpha 

%%
% No alpha channel is present, so |alpha| is empty.   

