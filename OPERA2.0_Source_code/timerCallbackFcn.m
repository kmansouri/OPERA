% This is the timer function to update the figure. 
% Save as timerCallbackFcn.m
function timerCallbackFcn(hTimer, eventData, hfig, I, map)
    figure(hfig);

    % i is the frame index
    persistent i;    
    if isempty(i), i=1; end

    % display the i-th frame in the GIF
    imshow(I(:,:,i),map);    

    % increment frame index i
    i=i+1;
    numframes=size(I,4);
    if (i>numframes), i=1; end
end