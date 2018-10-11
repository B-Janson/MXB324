function [] = CREATE_VIDEO( videofile, frames, framerate )

videofile.FrameRate = framerate;
open(videofile);

for i = 1:length(frames)
    writeVideo(videofile,frames(i));
end

close(videofile);

end