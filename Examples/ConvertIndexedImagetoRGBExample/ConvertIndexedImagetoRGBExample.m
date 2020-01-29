%% Convert Indexed Image to RGB  
% Read an image and convert it to an RGB image.   

% Copyright 2015 The MathWorks, Inc.


%% 
% Read the first image in the sample indexed image file, |corn.tif|. 
[X,map] = imread('corn.tif'); 

%%
% |X| is a 415-by-312 array of type |uint8|.  

%% 
% Verify that the colormap, |map|, is not empty, and convert the data in
% |X| to RGB. 
if ~isempty(map)
    Im = ind2rgb(X,map);
end  

%% 
% View the size and class of |X|. 
whos Im 

%%
% |X| is now a 415-by-312-by-3 array of type |double|.   

