function [FVCD] = update_brain_roi_colors_v2(roi_vals,roi_to_v_mat,input_struct,use_rwb)

    %roi_vals: nrois X 1 vector of values to be plotted for each ROI
    %roi_to_v_mat: vertices X roi matrix that transforms roi values into vertex values
    %
    %radius: a number that is used to determine the spread of the gaussian
    %function that is used to compute the contribution of each roi to each
    %vertex
    %
    %input_struct: to set the colorbar, set input_struct.clim = [cmin cmax]
    
    
    
    if ~exist('input_struct','var')
        input_struct = struct();
    end
    if isfield(input_struct,'clim')
        maxval = input_struct.clim(2);
        minval = input_struct.clim(1);
    else
        maxval = nanmax(roi_vals);
        minval = nanmin(roi_vals);
    end
    if ~exist('use_rwb','var')
        use_rwb = false;
    end
    
    if use_rwb
        load red_white_blue.mat
        ncolors = size(cmRedWhtBlu,1);
        cmap = cmRedWhtBlu;
    else
        ncolors = 1000;
        cmap = jet(ncolors);
    end
    
    %figure
    %[hs1 FVCD lights] = visualizeRawBrain(V,F);
    
    vertex_data = roi_to_v_mat*roi_vals(~isnan(roi_vals));%vertices X 1
    
    %keyboard
    num_v = numel(vertex_data);
    FVCD = repmat([.5 .5 .5],num_v,1);
    c_ind = round((vertex_data-minval)/(maxval-minval)*ncolors);
    non_nan_inds = ~isnan(c_ind);
    c_ind(c_ind<1)=1;
    c_ind(c_ind>ncolors) = ncolors;
    FVCD(non_nan_inds,:) = cmap(c_ind(non_nan_inds),:);
    
end