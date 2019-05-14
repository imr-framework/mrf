function Make_Movie(data)

dim1 = size(data);
v = VideoWriter('TMRF_v1.avi','Grayscale AVI');
colormap(gray(256));
v.FrameRate =10;
open(v);
k=1;


for t = 1:dim1(3)
    figure(1); 
    S = abs(data(:,:,t));
    imagesc(S(:,:)); 
    caxis([0 0.5.*max(S(:))]); colormap(gray);axis off;
    
  
    title(['Frame # ',num2str(t)]);
    pause(0.1);
    drawnow;
    S = S./max(abs(S(:)));
    writeVideo(v,S);
%     F(k) = getframe;
    k = k+1;
end
close (v);
end