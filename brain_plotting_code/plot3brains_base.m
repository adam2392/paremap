function [plots,brains] = plot3brains_base(i)
    close all
    figure(i)
    load brain.mat

    % gray color matrix
    FVCD = repmat([.5 .5 .5], size(V,1),1);
    views = [90 0; -180 -90; -90 0];

    numViews = size(views,1);

    brainCount = 0;
    plots = nan(3,1);
    brains = nan(3,1);
    set(gcf, 'Position', [200 200 1200 850])
    bot = .3;
    
    if ~exist('input_struct','var')
        input_struct = struct();
    end
    
    for v = 1:numViews,

        if v == 1
           pos = [-.065 .3 .5 .5];  % left bottom width height
        elseif v == 2
           pos = [.281 .33 .44 .44]; % left bottom width height   
        else
           pos = [.565 .3 .5 .5];  % left bottom width height   
        end
        
        % plot brains
        brainCount = brainCount + 1;
        
        plots(brainCount) = axes('Parent',gcf,'Position',pos);
        axis equal;
        % base brain
        brains(brainCount) = patch('faces',F,'vertices',V,'edgecolor','none','FaceColor','interp','FaceVertexCData',FVCD);
        
        set(gcf,'Renderer','OpenGL')
        set(gca,'visible','off');
        set(gca, 'cLim', [-3 3]);

        set(gcf,'Renderer','OpenGL')
        
        setBrainProps(brains(brainCount)); lighting phong
        view([0 0]);        camlight infinite;
        view([180 0]);      camlight infinite;
        view([-180 -90]);   camlight infinite;
        view([90 0]);       camlight infinite;
        view(views(v,:))
    end
end