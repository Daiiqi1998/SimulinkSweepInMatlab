classdef ModeSim < int32
% 0 - 不扫频
% 1 - 每组频点对应一次仿真（最准确）
% 2 - 所有频点合为一次长仿真（不推荐用于测量闭环振荡）
    enumeration  
        Signal(0)  
        Separate(1)  
        Together(2)  
    end  
end 