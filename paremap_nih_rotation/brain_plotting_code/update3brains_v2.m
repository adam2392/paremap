% FUNCTION: UPDATE3BRAINS
% Description: To update all 3 brain plots of front, side, top view
% correspondingly with color coded regions
% Input:
% - pat_s = patient ID (e.g. NIH039)
% - bp_flag = 0, or 1 for monopolar, or bipolar respectively
% Output:
% - els = struct containing metadata about each electrode
function h = update3brains_v2(brains,roi_vals,input_struct,title_str,cbar_title,use_rwb,h)
    try 
        delete(h);
    catch ploops
    end
    load red_white_blue.mat
    load ROI.mat
    load brain.mat
    
    ROItoplot = ROI(~isnan(roi_vals),:);
    roi_to_v_mat = create_roi_to_v_matrix(ROItoplot,V,12.5);

    views = [90 0; -180 -90; -90 0];

    numViews = size(views,1);

    brainCount = 0;
    bot = .3;
    
    if ~exist('input_struct','var')
        input_struct = struct();
    end
    if ~exist('use_rwb','var')
        use_rwb = false;
    end
    if isfield(input_struct,'clim')
        maxval = input_struct.clim(2);
        minval = input_struct.clim(1);
    else
        maxval = nanmax(roi_vals);
        minval = nanmin(roi_vals);
    end
    if isnan(minval) || isnan(maxval)
        input_struct.clim = [0 1];
    else
        input_struct.clim = [minval maxval];
    end
    
    FVCD = update_brain_roi_colors_v2(roi_vals,roi_to_v_mat,input_struct,use_rwb);

    for v = 1:numViews

        % plot brains
        brainCount = brainCount + 1;

        set(brains(brainCount),'FaceVertexCData',FVCD);

        if v == 2
             if use_rwb
                colormap(cmRedWhtBlu)
             else
                 colormap(jet)
             end
             h = colorbar('peer', gca, 'South');
             set(h, 'Position', [.3 bot-.075 .4 .03])
             set(h,'xaxislocation','bottom');
             set(gca, 'fontsize', 20);
             if exist('cbar_title','var')
                 set(get(h,'xlabel'),'string',cbar_title,'fontsize',20);
             end
             
             set(gca, 'cLim', input_struct.clim);
            
        end

        set(gcf,'Renderer','OpenGL')
        
    end
    if exist('title_str','var')
        s = suptitle(title_str);
        set(s,'Position',[.5 -.2 0]);
        set(s,'FontSize',30);
    end
end