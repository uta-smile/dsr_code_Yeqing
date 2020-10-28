classdef CascadeObjectDetector
    % overload subsref
    methods
        function ret = subsref(obj, I)
            ret = [0, 0, size(I.subs{1})];
        end
    end
end