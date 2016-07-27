% FUNCTION: ELECTOROI
% Description: Creates an ROI matrix from all the electrodes
%
% Input:
% - ROI = A matrix of ROIs x 3 with <x,y,z> coords of each ROI
% - elecs: A matrix of elecs x 3 with <x,y,z> coords of each electrode
% - radius: the radius of each ROI; all electrodes within this radius are in
% that ROI
% Output:
% - elecToROI = ROIs x elecs
function elecToROI = create_el_to_roi_matrix(ROI,elecs,radius)    
    dists = pdist2(ROI,elecs);%rois by elecs, distances between every elec and every roi
    valid_dists = dists < radius;%rois by elec, 1 if distance is less than some radius
    elecToROI = bsxfun(@rdivide,valid_dists,sum(valid_dists,2));%rois by elec, when mutliplied with elec by 1, will give values by ROI
end
