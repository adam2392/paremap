function h = plot_electrode_on_brains(plots,el,radius,color)
    for i = 1:3
        axes(plots(i));
        hold on
        h{i} = plot_electrode_v2(el,radius,color);
    end
end