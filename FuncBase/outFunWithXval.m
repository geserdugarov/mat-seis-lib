function stop = outFunWithXval(x, ~, state)
    stop = false;
    if strcmp(state,'iter')
        if size(x,1)==1
            disp(['                                                                              x   =   ',num2str(x)]);
        else
            disp(['                                                                              x   =   ',num2str(x')]);
        end
    end
end 