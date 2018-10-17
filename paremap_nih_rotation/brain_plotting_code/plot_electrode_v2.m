function h = plot_electrode_v2(el,radius,color)
    center = [el.x el.y el.z];
    %plotCircle3D(center,[0 0 0],radius,color);
    h(1) = plotCircle3D(center,[0 0 1],radius,color);
    h(2) = plotCircle3D(center,[0 1 0],radius,color);
    h(3) = plotCircle3D(center,[1 0 0],radius,color);
end