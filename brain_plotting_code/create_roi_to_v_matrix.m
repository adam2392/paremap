function roiToVert = create_roi_to_v_matrix(ROI,V,radius)    
    dists = pdist2(V,ROI);%vertices by rois, distances between every roi and every vertex
    valid_dists = dists < radius;%vertices by rois, 1 if distance is less than some radius
    weights = gauss(dists,(radius/3)^2).*valid_dists;%vertices by rois, weights to be applied to each vertex
    roiToVert = bsxfun(@rdivide,weights,sum(weights,2));
end

