function remove_electrode(h)
    for i = 1:3
        delete(h{i}(1));
        delete(h{i}(2));
        delete(h{i}(3));
    end
end