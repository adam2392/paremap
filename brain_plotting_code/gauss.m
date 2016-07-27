function x = gauss(dists,variance)
    x = exp(-(dists.^2)./(2*variance));
end