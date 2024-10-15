%% simulink 扫频程序 支持不同种类端口
% 
% 绪论：
%   当问题扩展到多端口时，需要测量的不只是单端口上的阻抗关系，可能还有单相电压对三相电压的传递关系。
%   所以这个程序在多模式扫频的基础上，允许了注入与测量的分开设置。
%   理论上也支持变步长仿真。
% 
% 使用说明：
%   只支持同时测量一组输入输出关系。请注意电流方向。
%   建议将 simulink 模型存放在与该程序同一路径。
%   它从 simulink 中获取数据依托于 slx 模型中的 to workspace "freq_swp" 模块，该模块的输入需要用户根据测量模式 INJ_MODE, MES_MODE 自己构造：
%       它需要有7组信号： inA, inB, inC, outA, outB, outC, theta.
%       如果需要测单相，那么请接到A相，BC用零占位。
%       theta 需要为弧度制。该程序仅在涉及dq模式下需要该信号，但即便你不需要它，也请用0占位，让输出中包含该维度.
%       （如果有界面，那么我会想用一个选项把 theta 合并至 inA 或者 outA. 可惜在这里写起来太麻烦了。）
%   如果文字说明不够形象，请看 example.slx .
%
% 运行所需支持文件：
%   ModeSim 和 ModeTest 在对应名字的 .m 文件里有定义。是为了方便维护代码，定义的枚举类。
%   myBodeCal 和 myBodeFix 也在对应名字的 .m 文件里有定义。是画图时候用到的方便函数。
% 
% 变量重名警告：
%   runT (1x1 double) 该次仿真模型总运行时间
%   injA (MxN double) 小信号扰动源的幅值，N为时域关键帧数，M为同时的扫频点数
%   injF (MxN double) 小信号扰动源的频率
%   injT (Mx1 double) 小信号扰动源开始作用的时间帧
%   injM (Mx1 int32) 注入模式
%   请不要在模型中重新定义（比如InitFcn）这些参数，或是将它们应用在扫频模块之外的模型中。另外，请记得把模型总运行时间设为runT。

%% 扫频设置
% 如果你是用户，请你只改动这部分的变量，没有特殊需求就不要改其它部分了。

model_name = "example";    % 模型名
save_name = "Swp_A1";    % 扫频结果保存名。若缺省，则按当前时间保存。后面保存的时候会自动加上模式注释。
% 如果还需要调用别的路径下的自定义函数，请在这里 addpath

SWP_MODE = ModeSim.Separate;   % 仿真模式：
% 下述“注入端口”不必须是仿真中真正注入谐波的端口。但在三相双轴测试模式下，需要保证在谐波源单轴注入时，“注入端口”的谐波没有混杂。
INJ_MODE = ModeTest.P3_pn;     % 注入端口测量模式
MES_MODE = ModeTest.P3_pn;     % 测量端口测量模式


MES_COUP = 0;   % 测量端口相对注入端口的基准频率耦合情况。
% 即，测量多端口频率耦合时，测量端口使用的 F_FREQS' 是否需要做特殊位移。 
% -2, -1, 0, 1, 2 -- "0"表示同样频率。在此基础上，每个数表示测量端口额外加减的基波频率个数。即， F_FREQS' = F_FREQS + MES_COUP * F_FUND。
% （备用功能）98, 99, 100, 101, 102 -- "100"表示测量端口是注入端口关于基波的对称。在此基础上，加减基波个数。即， F_FREQS' = (2*F_FUND-F_FREQS) + (MES_COUP-100) * F_FUND。
% 注：正负序测量模式下，由于对每组注入都会进行两次测量，所以将本征频率（而非耦合频率）当作 MES_COUP 所描述的端口频率。也就是，如果你想用 ModeTest.P3_pn 选项测量同一端口的正负序阻抗矩阵，那么请设置 MES_COUP = 0， 而不是100。

AUTO_MODE = true;  % AUTO模式。若启动该选项，则认为注入、测量的端口性质相同，自动设置 MES_MODE = INJ_MODE, MES_COUP = 0。

INJ_NUM = 5;    % 允许同时注入的最大频点数。因为仿真里的谐波是在遍历函数里生成的，所以理论上不用设限。

% round( logspace( log10(最小频率),log10(最大频率),取点数量 )/频率分辨率 )*频率分辨率
F_FREQS = unique([...
%    round( logspace(log10(5),   log10(40), 10)/1  )*1, ...
%    round( logspace(log10(40),   log10(100), 15)/5  )*5, ...
%    round( logspace(log10(100), log10(5e3), 30)/25  )*25])';
round( logspace(log10(5), log10(5e3), 5)/10  )*10])';
% 该频率是注入端口的频率。注意请把它设置为竖的向量。
% 在dq模式下，它是在dq系下的等效注入频率。
% 如果你想注入负序，请直接用负频率代替。若它同时包含正频率与负频率，下面画图时会分出一个正半图，一个负半图。
% 在正负序，AUTO模式下，它会取该频率的绝对值作为正序注入集，并自动生成负序。

T_SIM_MIN = 1e-6;   % 仿真中的时间最小值（该扫频程序允许离散与连续仿真）

F_FUND = 50;    % 基波频率
F_TICK_MAX = 10;    % 频率最大分辨率。它防止高频测试中忽略基波和纹波的影响。严格来说，请把它设成基波频率的因数。
injA_ref = 1;   % 注入小信号幅值

T_STEADY = 0.2;   % 认为从此开始，仿真模型处于稳态
T_FADE = 0.1;  % ModeSim.Together情况下使用。同一仿真中，切换小信号时的等待时间。被测系统增益越高，建议等待越长（暂时没做实时自适应）

%% 对仿真设置的处理

if MES_MODE == ModeTest.P3_dq
    INJ_NUM = max(INJ_NUM,2);
end

% 去零
F_FREQS(F_FREQS == 0) = [];
f_freqs_pos = sort(F_FREQS(F_FREQS>0));
f_freqs_neg = -sort(-(F_FREQS(F_FREQS < 0)));
F_FREQS = [f_freqs_pos;f_freqs_neg];    % 排版：先是从小到大排序的正数部分，再是绝对值从小到大排序的负数部分

% 频率分辨率
f_dec = F_FREQS(F_FREQS ~= fix(F_FREQS));
if ~isempty(f_dec)
    % 有小数，向小于1方向搜索
    temp_10 = 0;
    while any(mod(f_dec,1) ~= 0)
        f_dec = f_dec*10;
        temp_10 = temp_10 - 1;
    end
    f_tick_min = exp(10, temp_10);
    clear temp_10
else
    % 全整数，向大于1方向搜索
    f_tick_min = F_FREQS(1);
    for i = 2:length(F_FREQS)
        f_tick_min = gcd(f_tick_min, F_FREQS(i));
    end
    isOnly2And5 = (f_tick_min == 1);
    if ~isOnly2And5
        f_tick_copy = f_tick_min;
        while mod(f_tick_copy, 2) == 0
            f_tick_copy = f_tick_copy / 2;
        end
        while mod(f_tick_copy, 5) == 0
            f_tick_copy = f_tick_copy / 5;
        end
        isOnly2And5 = (f_tick_copy == 1);
    end

    while ~isOnly2And5
        f_tick_min = f_tick_min - 1;

        isOnly2And5 = (f_tick_min == 1);
        if ~isOnly2And5
            f_tick_copy = f_tick_min;
            while mod(f_tick_copy, 2) == 0
                f_tick_copy = f_tick_copy / 2;
            end
            while mod(f_tick_copy, 5) == 0
                f_tick_copy = f_tick_copy / 5;
            end
            isOnly2And5 = (f_tick_copy == 1);
        end
        
    end
    clear isOnly2And5 f_tick_copy
end
clear f_dec

% 如果选了自动，则接管一些设置
if AUTO_MODE
    MES_MODE = INJ_MODE;
    MES_COUP = 0;
    F_FREQS = sort(unique(abs(F_FREQS)));   % 去负
    if INJ_MODE == ModeTest.P3_pn
        for i=1:length(F_FREQS)
            if F_FREQS(i)<=2*F_FUND
                F_FREQS = [F_FREQS;-F_FREQS(i)];
            else
                F_FREQS = [F_FREQS;2*F_FUND - F_FREQS(i)];
            end
        end
        F_FREQS(F_FREQS == 2*F_FUND) = [];
    else
        if INJ_MODE == ModeTest.P3_dq
            F_FREQS(F_FREQS == F_FUND) = [];
        end
    end
end
% figure(1); semilogx(F_FREQS,ones(1,length(F_FREQS)),"*")   % 画个图验证扫频点排布，可以注释掉
N_FREQS = length(F_FREQS);

% 非自动的情况下，有些设置也可以接管。比如，若非三相正负序，则负频率无效（三相负频率等于负序，两相负频率等于共轭，单相负频率等于互补）。
if INJ_MODE ~= ModeTest.P3_pn
    F_FREQS = unique(abs(F_FREQS));
end

% 由于 F_FREQS 可能发生了变化，所以重新排版（排版信息在后续换算结果时会使用到）
f_freqs_pos = sort(F_FREQS(F_FREQS>0));
f_freqs_neg = -sort(-(F_FREQS(F_FREQS < 0)));
F_FREQS = [f_freqs_pos;f_freqs_neg];

% 生成实际的测量频点（中心频点）
if MES_COUP < 50
    F_MEASU = F_FREQS + MES_COUP * F_FUND;
else
    F_MEASU = (2*F_FUND-F_FREQS) + (MES_COUP-100) * F_FUND;
end

F_MEASU_COUP = 2*F_FUND-F_MEASU;    % 仅为 P3_pn 模式的测量准备

% 去掉无效频点（F_FREQS, F_MEASU, F_MEASU_COUP 中任何一个可能为0的项）
F_MEASU_COUP(F_FREQS == 0) = [];
F_MEASU(F_FREQS == 0) = [];
F_FREQS(F_FREQS == 0) = [];
F_MEASU_COUP(F_MEASU == 0) = [];
F_FREQS(F_MEASU == 0) = [];
F_MEASU(F_MEASU == 0) = [];
if INJ_MODE == ModeTest.P3_pn
F_FREQS(F_MEASU_COUP == 0) = [];
F_MEASU(F_MEASU_COUP == 0) = [];
F_MEASU_COUP(F_MEASU_COUP == 0) = [];
end
if ~(numel(F_FREQS) == numel(F_MEASU) && numel(F_FREQS) == numel(F_MEASU_COUP))
    disp('频率去重部分出错，请联系调试')
    return
end

k_inj = F_FREQS/F_FUND;
k_mes = F_MEASU/F_FUND;

% 生成最大的测量频点
f_temp = [F_MEASU; F_MEASU_COUP];
f_mes_abs_max = max(max(F_FREQS), -min(F_FREQS));
if INJ_MODE == ModeTest.P3_dq
    f_mes_abs_max = max(f_mes_abs_max, max(F_FREQS) + F_FUND);
elseif INJ_MODE == ModeTest.P3_pn
    f_mes_abs_max = max(f_mes_abs_max, 2*F_FUND - min(F_FREQS));
end
f_mes_abs_max = max([f_mes_abs_max, max(f_temp), -min(f_temp)]);
if INJ_MODE == ModeTest.P3_dq
    f_mes_abs_max = max(f_mes_abs_max, max(f_temp) + F_FUND);
elseif INJ_MODE == ModeTest.P3_pn
    f_mes_abs_max = max(f_mes_abs_max, 2*F_FUND - min(f_temp));
end

f_mes_max = 1;
while (f_mes_max*10 <= f_mes_abs_max)
    f_mes_max = f_mes_max*10;
end
if (f_mes_max*2 >= f_mes_abs_max)
    f_mes_max = f_mes_max*2;
elseif (f_mes_max*5 >= f_mes_abs_max)
    f_mes_max = f_mes_max*5;
else
    f_mes_max = f_mes_max*10;
end

T_TICK = min(T_SIM_MIN, 0.1/f_mes_max); % 这是时间的最小分辨率。出于测量需求，它可以比仿真分辨率更小。
if T_SIM_MIN > 1/f_mes_max
    warning("仿真时长大于测量分辨率，高频结果不可靠")
elseif T_SIM_MIN > 0.1/f_mes_max
    warning("仿真时长小于测量分辨率，但大于1/10测量分辨率，高频结果可能有偏差");
end
clear f_mes_abs_max f_mes_max f_temp

% 频点分组
if MES_MODE == ModeTest.P3_dq
    inj_num = floor(INJ_NUM/2);
else
    inj_num = INJ_NUM;
end
f_freqs_group = zeros(1,inj_num);
f_measu_group = zeros(1,inj_num);
f_measu_coup_group = zeros(1,inj_num);
index_group = zeros(1,inj_num);
bool_freq = true(1, N_FREQS);
i = 1;
while any(bool_freq)
    f_freqs_group(i,:) = zeros(1,inj_num);  % 占位新行
    f_measu_group(i,:) = zeros(1,inj_num);
    f_measu_coup_group(i,:) = zeros(1,inj_num); % 占位新行
    temp_index = find(bool_freq, 1);
    temp_row = [];
    j = 1;
    for k = temp_index : N_FREQS
        if bool_freq(k)&&(isempty(f_freqs_group(i,:)) || ~any(ismember(temp_row, F_FREQS(k))))
            f_freqs_group(i,j) = F_FREQS(k);
            f_measu_group(i,j) = F_MEASU(k);
            f_measu_coup_group(i,j) = F_MEASU_COUP(k);
            index_group(i,j) = k;
            j = j+1;
            bool_freq(k) = false;
            temp_row = [temp_row, F_FREQS(k)];
            if ~(INJ_MODE == ModeTest.P1 && MES_MODE == ModeTest.P1)
                temp_row = [temp_row, -F_FREQS(k), ...  % 同一频率可能的正负半轴测量结果
                            F_FREQS(k) + F_FUND, F_FREQS(k) - F_FUND, ...   % 传递环节中可能的谐波耦合
                            2*F_FUND - F_FREQS(k)];     % 传递环节中可能的正负序对称
            end
            if MES_COUP < 50
                temp_row = [temp_row, F_FREQS(k) + MES_COUP * F_FUND];      % MES_COUP 决定的特殊耦合项
            else
                temp_row = [temp_row, (2*F_FUND-F_FREQS(k)) + MES_COUP * F_FUND];
            end
            temp_row = unique(temp_row);
        end
        if j > inj_num
            break
        end
    end
    i = i+1;
end
[f_group_m, f_group_n] = size(f_freqs_group);

% 分组后频点的时间分辨率计算
f_ticks_group = zeros(f_group_m,1);
for i = 1:f_group_m
    f_tick = f_tick_min;
    f_freqs_row = f_freqs_group(i, :);
    for j = 1:length(f_freqs_row)
        while all( mod(f_freqs_row/(f_tick*2),1)==0 )
            f_tick = f_tick * 2;
        end
        while all( mod(f_freqs_row/(f_tick*5),1)==0 )
            f_tick = f_tick * 5;
        end
    end
    f_tick = min(f_tick, F_TICK_MAX);
    f_ticks_group(i) = f_tick;
end

% 占缓存
if INJ_MODE == ModeTest.P3_ab || INJ_MODE == ModeTest.P3_dq
    % 需要同频不同轴单独注入的情况下，不会有额外的频点，需要人工控制扫两轮
    M_FREQS = 2*f_group_m;  % 需要单独注入的次数
    Comp_i = zeros(N_FREQS, 2);     % 行内存储不同轴注入时的结果
    Comp_I = zeros(N_FREQS, 2);     % 行内存储稳态值在不同轴的分量
    Comp_o_self = zeros(N_FREQS, 2);    
    Comp_o_coup = zeros(N_FREQS, 2);    % 如果待测量目标不是三相（二维）的，那么不会用到它。
else
    M_FREQS = f_group_m;
    Comp_i = zeros(N_FREQS, 1);
    Comp_I = zeros(N_FREQS, 1);
    Comp_o_self = zeros(N_FREQS, 1);
    Comp_o_coup = zeros(N_FREQS, 1);
end

Comp_O_ax1 = zeros(N_FREQS, 1);
Comp_O_ax2 = zeros(N_FREQS, 1);

%{
% 给程序维护人员看的：截止到这里，我们有了这些参数：
% 扫频模式 SWP_MODE = ModeSim.Signal / Separate / Together, 注入端、输出端测量模式 INJ_MODE, MES_MODE = ModeTest.P1 / P3_ab / P3_dq / P3_pn 
% 两个端口间频率耦合情况 MES_COUP = -2/-1/0/1/2/100, 同时注入频点数 INJ_NUM = 5, 注入幅值 injA_ref
% 注入频率 F_FREQS = [f_freqs_pos; f_freqs_neg], 其个数 N_FREQS 与次数 k_inj, 测量频率 F_MEASU 及其次数 k_mes, 基波频率 F_FUND
% 分好组的注入频率 f_freqs_group 与输出频率 f_measu_group, f_measu_coup_group, 伴生的 index_group, 它所对应的实际注入次数 M_FREQS
% 时间步长 T_TICK, 仿真等待时间 T_STEADY, 小信号等待时间 T_FADE
% 频率最小分辨率 f_tick_min, 最大分辨率 F_TICK_MAX, 与注入频率分组对应的分辨率 f_ticks_group
%}


%% 第一次运行仿真：获取稳态工作点矢量
f_tick = f_tick_min;
T_SCALE = 1/f_tick;     % 根据取频点的分辨率计算需要取的时间段长度
% T_SCALE = round(1/K_PER)/F_FUND;

% Signal 模式，构造设置并仿真
simIn = Simulink.SimulationInput(model_name);
simIn = setVariable(simIn, "runT", T_STEADY+T_SCALE);
simIn = setVariable(simIn, "injA", [0,]);
simIn = setVariable(simIn, "injF", [50,]);
simIn = setVariable(simIn, "injT", [0,]);
simIn = setVariable(simIn, "injM", int32(ModeTest.signal));
simOut = sim(simIn);
[result_m, result_n] = size(simOut.freq_swp.signals.values);      % 用于验证维度是否与模式设置匹配
if result_n < 7
    disp("输出通道数不足");
    return
end

% 不管仿真是定步长、变步长，对采样线性插值就行
t = T_STEADY : T_TICK : simOut.freq_swp.time(end);      N_t = length(t);    % 取稳态时间段


switch INJ_MODE
    case ModeTest.P1
        % 获取时域测量，采样插值以模拟定点示波器
        x = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
        y = fft(x)*2/N_t;
        
        % 提取响应
        for i=1:N_FREQS
            if F_FREQS(i) >= 0
                Comp_I(i,1) = y(round(1 + F_FREQS(i)/f_tick));
            else
                Comp_I(i,1) = y(round(end+1 + F_FREQS(i)/f_tick));
            end
        end

    case ModeTest.P3_ab
        x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
        x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
        x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
        x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
        x_beta = sqrt(1/2) * (x_b - x_c);
        y_alpha = fft(x_alpha)*2/N_t;
        y_beta = fft(x_beta)*2/N_t;

        for i=1:N_FREQS
            if F_FREQS(i) >= 0
                Comp_I(i,1) = y_alpha(round(1 + F_FREQS(i)/f_tick));
                Comp_I(i,2) = y_beta(round(1 + F_FREQS(i)/f_tick));
            else
                Comp_I(i,1) = y_alpha(round(end+1 - F_FREQS(i)/f_tick));
                Comp_I(i,2) = y_beta(round(end+1 - F_FREQS(i)/f_tick));
            end
        end

    case ModeTest.P3_dq
        x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
        x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
        x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
        theta = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,7), t);
        cos_t = cos(theta);
        sin_t = sin(theta);
        x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
        x_beta = sqrt(1/2) * (x_b - x_c);
        x_d = cos_t .* x_alpha + sin_t .* x_beta;
        x_q = -sin_t .* x_alpha + cos_t .* x_beta;
        y_d = fft(x_d)*2/N_t;
        y_q = fft(x_q)*2/N_t;
        
        for i=1:N_FREQS
            if F_FREQS(i) >= 0
                Comp_I(i,1) = y_d(round(1 + F_FREQS(i)/f_tick));
                Comp_I(i,2) = y_q(round(1 + F_FREQS(i)/f_tick));
            else
                Comp_I(i,1) = y_d(round(end+1 - F_FREQS(i)/f_tick));
                Comp_I(i,2) = y_q(round(end+1 - F_FREQS(i)/f_tick));
            end
        end
        
    case ModeTest.P3_pn
        x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
        x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
        x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
        y_a = fft(x_a)*2/N_t;
        y_b = fft(x_b)*2/N_t;
        y_c = fft(x_c)*2/N_t;
        a = exp(1j*2*pi/3); % 相序分解需要的旋转常量
        y = (y_a + a*y_b + a*a*y_c)/3;
        
        for i=1:N_FREQS
            if F_FREQS(i) >= 0
                Comp_I(i,1) = y(round(1 + F_FREQS(i)/f_tick));
            else
                Comp_I(i,1) = y(round(end+1 + F_FREQS(i)/f_tick));
            end
        end

    otherwise
        disp("注入端口模式违法");
        return
end

switch MES_MODE
    case ModeTest.P1
        % 获取时域测量，采样插值以模拟定点示波器
        x = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
        y = fft(x)*2/N_t;
        
        % 提取响应
        for i=1:length(F_MEASU)
            if F_MEASU(i) >= 0
                Comp_O_ax1(i,1) = y(round(1 + F_MEASU(i)/f_tick));
            else
                Comp_O_ax1(i,1) = y(round(end+1 + F_MEASU(i)/f_tick));
            end
        end

    case ModeTest.P3_ab
        x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
        x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
        x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
        x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
        x_beta = sqrt(1/2) * (x_b - x_c);
        y_alpha = fft(x_alpha)*2/N_t;
        y_beta = fft(x_beta)*2/N_t;

        for i=1:length(F_MEASU)
            if F_MEASU(i) >= 0
                Comp_O_ax1(i,1) = y_alpha(round(1 + F_MEASU(i)/f_tick));
                Comp_O_ax2(i,1) = y_beta(round(1 + F_MEASU(i)/f_tick));
            else
                Comp_O_ax1(i,1) = y_alpha(round(end+1 - F_MEASU(i)/f_tick));
                Comp_O_ax2(i,1) = y_beta(round(end+1 - F_MEASU(i)/f_tick));
            end
        end

    case ModeTest.P3_dq
        x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
        x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
        x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
        theta = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,7), t);
        cos_t = cos(theta);
        sin_t = sin(theta);
        x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
        x_beta = sqrt(1/2) * (x_b - x_c);
        x_d = cos_t .* x_alpha + sin_t .* x_beta;
        x_q = -sin_t .* x_alpha + cos_t .* x_beta;
        y_d = fft(x_d)*2/N_t;
        y_q = fft(x_q)*2/N_t;
        
        for i=1:length(F_MEASU)
            if F_MEASU(i) >= 0
                Comp_O_ax1(i,1) = y_d(round(1 + F_MEASU(i)/f_tick));
                Comp_O_ax2(i,1) = y_q(round(1 + F_MEASU(i)/f_tick));
            else
                Comp_O_ax1(i,1) = y_d(round(end+1 - F_MEASU(i)/f_tick));
                Comp_O_ax2(i,1) = y_q(round(end+1 - F_MEASU(i)/f_tick));
            end
        end
        
    case ModeTest.P3_pn
        x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
        x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
        x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
        y_a = fft(x_a)*2/N_t;
        y_b = fft(x_b)*2/N_t;
        y_c = fft(x_c)*2/N_t;
        a = exp(1j*2*pi/3); % 相序分解需要的旋转常量
        y = (y_a + a*y_b + a*a*y_c)/3;
        
        for i=1:length(F_MEASU)
            if F_MEASU(i) >= 0
                Comp_O_ax1(i,1) = y(round(1 + F_MEASU(i)/f_tick));
            else
                Comp_O_ax1(i,1) = y(round(end+1 + F_MEASU(i)/f_tick));
            end
            if F_MEASU_COUP(i) >= 0
                Comp_O_ax2(i,1) = y(round(1 + F_MEASU_COUP(i)/f_tick));
            else
                Comp_O_ax2(i,1) = y(round(end+1 + F_MEASU_COUP(i)/f_tick));
            end
        end

    otherwise
        disp("输出端口模式违法");
        return
end



%% 遍历运行仿真：测量小信号响应

switch SWP_MODE
    case ModeSim.Signal
        disp("已经结束了。")
        return
    case ModeSim.Separate
    
        % Separate 模式，循环扫频
        
        for ix = 1:M_FREQS

            i = mod(ix-1, f_group_m)+1;    % 将ix折算到1~f_group_m间
            f_freqs_row = f_freqs_group(i, :);      % 取当前注入行
            f_measu_row = f_measu_group(i, :);
            f_measu_coup_row = f_measu_coup_group(i, :);
            index_row = index_group(i, :);

            % 计算当前注入行对应的分辨率
            f_tick = f_ticks_group(i);
            T_SCALE = 1/f_tick;

            % 决定当前注入行的模式
            switch INJ_MODE
                case ModeTest.P1
                    injM_ref = ModeTest.P1;
                case ModeTest.P3_ab
                    if ix <= f_group_m
                        injM_ref = ModeTest.P3_a;
                    else
                        injM_ref = ModeTest.P3_b;
                    end
                case ModeTest.P3_dq
                    if ix <= f_group_m
                        injM_ref = ModeTest.P3_d;
                    else
                        injM_ref = ModeTest.P3_q;
                    end
                case ModeTest.P3_pn
                    injM_ref = ModeTest.P3_pn;
                otherwise
                    injM_ref = ModeTest.signal;
            end

            % 构造设置并仿真
            simIn = Simulink.SimulationInput(model_name);
            simIn = setVariable(simIn, "runT", T_STEADY + T_SCALE);
            simIn = setVariable(simIn, "injA", injA_ref * ones(1,f_group_n));
            simIn = setVariable(simIn, "injF", f_freqs_row);
            simIn = setVariable(simIn, "injT", [0,]);
            
            simIn = setVariable(simIn, "injM", int32(injM_ref));
            simOut = sim(simIn);

            % 提取结果
            t = T_STEADY : T_TICK : simOut.freq_swp.time(end);      N_t = length(t);

            switch INJ_MODE
                case ModeTest.P1

                    x = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    y = fft(x)*2/N_t;
                    
                    % 提取响应
                    for j=1:f_group_n
                        if f_freqs_row(j) == 0
                            break
                        elseif f_freqs_row(j) > 0
                            Comp_i(index_row(j),1) = y(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        else
                            Comp_i(index_row(j),1) = y(round(end+1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        end
                    end

                case ModeTest.P3_ab

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    y_alpha = fft(x_alpha)*2/N_t;
                    y_beta = fft(x_beta)*2/N_t;

                    if ix <= f_group_m
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),1) = y_alpha(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            else
                                Comp_i(index_row(j),1) = y_alpha(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            end
                        end
                        
                    else
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),2) = y_beta(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            else
                                Comp_i(index_row(j),2) = y_beta(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            end
                        end
                        
                    end

                case ModeTest.P3_dq

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
                    theta = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,7), t);
                    cos_t = cos(theta);
                    sin_t = sin(theta);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    x_d = cos_t .* x_alpha + sin_t .* x_beta;
                    x_q = -sin_t .* x_alpha + cos_t .* x_beta;
                    y_d = fft(x_d)*2/N_t;
                    y_q = fft(x_q)*2/N_t;
                    
                    if ix <= f_group_m
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),1) = y_d(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            else
                                Comp_i(index_row(j),1) = y_d(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            end
                        end
                        
                    else
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),2) = y_q(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            else
                                Comp_i(index_row(j),2) = y_q(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            end
                        end
                        
                    end

                case ModeTest.P3_pn
                    
                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
                    y_a = fft(x_a)*2/N_t;
                    y_b = fft(x_b)*2/N_t;
                    y_c = fft(x_c)*2/N_t;
                    a = exp(1j*2*pi/3); % 相序分解需要的旋转常量
                    y = (y_a + a*y_b + a*a*y_c)/3;
                    
                    for j=1:f_group_n
                        if f_freqs_row(j) == 0
                                break
                        elseif f_freqs_row(j) > 0
                            Comp_i(index_row(j),1) = y(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        else
                            Comp_i(index_row(j),1) = y(round(end+1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        end
                    end

            end

            switch MES_MODE
                case ModeTest.P1

                    % 获取时域测量，采样插值以模拟定点示波器
                    x = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    y = fft(x)*2/N_t;
                    
                    % 提取响应
                    if ix <= f_group_m

                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                        end

                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                        end
                    
                    end

                case ModeTest.P3_ab

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    y_alpha = fft(x_alpha)*2/N_t;
                    y_beta = fft(x_beta)*2/N_t;

                    if ix <= f_group_m
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y_alpha(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_beta(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y_alpha(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_beta(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y_alpha(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_beta(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y_alpha(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_beta(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    end

                case ModeTest.P3_dq

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
                    theta = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,7), t);
                    cos_t = cos(theta);
                    sin_t = sin(theta);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    x_d = cos_t .* x_alpha + sin_t .* x_beta;
                    x_q = -sin_t .* x_alpha + cos_t .* x_beta;
                    y_d = fft(x_d)*2/N_t;
                    y_q = fft(x_q)*2/N_t;
                    
                    if ix <= f_group_m

                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y_d(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_q(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y_d(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_q(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end

                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y_d(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_q(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y_d(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_q(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    end

                case ModeTest.P3_pn
                    
                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
                    y_a = fft(x_a)*2/N_t;
                    y_b = fft(x_b)*2/N_t;
                    y_c = fft(x_c)*2/N_t;
                    a = exp(1j*2*pi/3); % 相序分解需要的旋转常量
                    y = (y_a + a*y_b + a*a*y_c)/3;
                    
                    if ix <= f_group_m

                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                            if f_measu_coup_row(j) == 0
                                break
                            elseif f_measu_coup_row(j) > 0
                                Comp_o_coup(index_row(j),1) = y(round(1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_coup(index_row(j),1) = y(round(end+1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end

                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                            if f_measu_coup_row(j) == 0
                                break
                            elseif f_measu_coup_row(j) > 0
                                Comp_o_coup(index_row(j),2) = y(round(1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_coup(index_row(j),2) = y(round(end+1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    end

            end
        
        end
        
    case ModeSim.Together
    
        % Together 模式，单次仿真

        % 构造设置
        t_temp = 0;
        runT_ref = T_STEADY - T_FADE;
        injT_ref = zeros(M_FREQS, 1);   % 开始注入的时间点
        mesT_ref = zeros(M_FREQS, 1);   % 开始测量的时间点。从这里开始测量一个 T_SCALE
        injM_ref = zeros(M_FREQS, 1);
        for ix = 1:M_FREQS

            % 每个频点的注入时长
            f_tick = f_ticks_group(mod(ix-1,f_group_m)+1);
            T_SCALE = 1/f_tick;
        
            % 每个频点的注入与测量时间点
            injT_ref(ix) = t_temp;
            if ix == 1
                mesT_ref(ix) = T_STEADY;
                t_temp = T_STEADY + T_SCALE;
            else
                mesT_ref(ix) = t_temp + T_FADE;
                t_temp = t_temp + T_FADE + T_SCALE;
            end
            
            runT_ref = runT_ref + T_FADE + T_SCALE;
        
            % 决定每个频点的注入模式
            switch INJ_MODE
                case ModeTest.P1
                    injM_ref(ix) = ModeTest.P1;
                case ModeTest.P3_ab
                    if ix <= f_group_m
                        injM_ref(ix) = ModeTest.P3_a;
                    else
                        injM_ref(ix) = ModeTest.P3_b;
                    end
                case ModeTest.P3_dq
                    if ix <= f_group_m
                        injM_ref(ix) = ModeTest.P3_d;
                    else
                        injM_ref(ix) = ModeTest.P3_q;
                    end
                case ModeTest.P3_pn
                    injM_ref(ix) = ModeTest.P3_pn;
                otherwise
                    injM_ref(ix) = ModeTest.signal;
            end
        end

        % 仿真
        simIn = Simulink.SimulationInput(model_name);
        simIn = setVariable(simIn, "runT", runT_ref);
        simIn = setVariable(simIn, "injA", injA_ref * ones(M_FREQS, f_group_n));
        if M_FREQS <= f_group_m
            simIn = setVariable(simIn, "injF", f_freqs_group);
        else
            simIn = setVariable(simIn, "injF", [f_freqs_group; f_freqs_group]);
        end
        simIn = setVariable(simIn, "injT", injT_ref);
        simIn = setVariable(simIn, "injM", injM_ref);
        simOut = sim(simIn);

        % 提取结果
        for ix = 1:M_FREQS

            i = mod(ix-1, f_group_m)+1;
            f_freqs_row = f_freqs_group(i, :);      % 取当前注入行
            f_measu_row = f_measu_group(i, :);
            f_measu_coup_row = f_measu_coup_group(i, :);
            index_row = index_group(i, :);
            
            f_tick = f_ticks_group(i);
            T_SCALE = 1/f_tick;
            
            t = mesT_ref(ix) : T_TICK : mesT_ref(ix) + T_SCALE;       N_t = length(t);
            
            switch INJ_MODE
                case ModeTest.P1

                    x = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    y = fft(x)*2/N_t;
                    
                    % 提取响应
                    for j=1:f_group_n
                        if f_freqs_row(j) == 0
                            break
                        elseif f_freqs_row(j) > 0
                            Comp_i(index_row(j),1) = y(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        else
                            Comp_i(index_row(j),1) = y(round(end+1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        end
                    end

                case ModeTest.P3_ab

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    y_alpha = fft(x_alpha)*2/N_t;
                    y_beta = fft(x_beta)*2/N_t;

                    if ix <= f_group_m
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),1) = y_alpha(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            else
                                Comp_i(index_row(j),1) = y_alpha(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            end
                        end
                        
                    else
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),2) = y_beta(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            else
                                Comp_i(index_row(j),2) = y_beta(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            end
                        end
                        
                    end

                case ModeTest.P3_dq

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
                    theta = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,end), t);
                    cos_t = cos(theta);
                    sin_t = sin(theta);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    x_d = cos_t .* x_alpha + sin_t .* x_beta;
                    x_q = -sin_t .* x_alpha + cos_t .* x_beta;
                    y_d = fft(x_d)*2/N_t;
                    y_q = fft(x_q)*2/N_t;
                    
                    if ix <= f_group_m
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),1) = y_d(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            else
                                Comp_i(index_row(j),1) = y_d(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                            end
                        end
                        
                    else
                    
                        for j=1:f_group_n
                            if f_freqs_row(j) == 0
                                break
                            elseif f_freqs_row(j) > 0
                                Comp_i(index_row(j),2) = y_q(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            else
                                Comp_i(index_row(j),2) = y_q(round(end+1 - f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),2);
                            end
                        end
                        
                    end

                case ModeTest.P3_pn
                    
                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,1), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,2), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,3), t);
                    y_a = fft(x_a)*2/N_t;
                    y_b = fft(x_b)*2/N_t;
                    y_c = fft(x_c)*2/N_t;
                    a = exp(1j*2*pi/3); % 相序分解需要的旋转常量
                    y = (y_a + a*y_b + a*a*y_c)/3;
                    
                    for j=1:f_group_n
                        if f_freqs_row(j) == 0
                                break
                        elseif f_freqs_row(j) > 0
                            Comp_i(index_row(j),1) = y(round(1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        else
                            Comp_i(index_row(j),1) = y(round(end+1 + f_freqs_row(j)/f_tick)) - Comp_I(index_row(j),1);
                        end
                    end

            end

            switch MES_MODE
                case ModeTest.P1

                    % 获取时域测量，采样插值以模拟定点示波器
                    x = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    y = fft(x)*2/N_t;
                    
                    % 提取响应
                    if ix <= f_group_m

                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                        end

                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                        end
                    
                    end

                case ModeTest.P3_ab

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    y_alpha = fft(x_alpha)*2/N_t;
                    y_beta = fft(x_beta)*2/N_t;

                    if ix <= f_group_m
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y_alpha(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_beta(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y_alpha(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_beta(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y_alpha(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_beta(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y_alpha(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_beta(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    end

                case ModeTest.P3_dq

                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
                    theta = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,end), t);
                    cos_t = cos(theta);
                    sin_t = sin(theta);
                    x_alpha = sqrt(2/3) * (x_a - 1/2 * x_b - 1/2 * x_c);
                    x_beta = sqrt(1/2) * (x_b - x_c);
                    x_d = cos_t .* x_alpha + sin_t .* x_beta;
                    x_q = -sin_t .* x_alpha + cos_t .* x_beta;
                    y_d = fft(x_d)*2/N_t;
                    y_q = fft(x_q)*2/N_t;
                    
                    if ix <= f_group_m

                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y_d(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_q(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y_d(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),1) = y_q(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end

                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y_d(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_q(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y_d(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                                Comp_o_coup(index_row(j),2) = y_q(round(end+1 - f_measu_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    end

                case ModeTest.P3_pn
                    
                    x_a = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,4), t);
                    x_b = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,5), t);
                    x_c = interp1(simOut.freq_swp.time, simOut.freq_swp.signals.values(:,6), t);
                    y_a = fft(x_a)*2/N_t;
                    y_b = fft(x_b)*2/N_t;
                    y_c = fft(x_c)*2/N_t;
                    a = exp(1j*2*pi/3); % 相序分解需要的旋转常量
                    y = (y_a + a*y_b + a*a*y_c)/3;
                    
                    if ix <= f_group_m

                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),1) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),1) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                            if f_measu_coup_row(j) == 0
                                break
                            elseif f_measu_coup_row(j) > 0
                                Comp_o_coup(index_row(j),1) = y(round(1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_coup(index_row(j),1) = y(round(end+1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                            
                            
                        end

                    else
                    
                        for j=1:f_group_n
                            if f_measu_row(j) == 0
                                break
                            elseif f_measu_row(j) > 0
                                Comp_o_self(index_row(j),2) = y(round(1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            else
                                Comp_o_self(index_row(j),2) = y(round(end+1 + f_measu_row(j)/f_tick)) - Comp_O_ax1(index_row(j),1);
                            end
                            if f_measu_coup_row(j) == 0
                                break
                            elseif f_measu_coup_row(j) > 0
                                Comp_o_coup(index_row(j),2) = y(round(1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            else
                                Comp_o_coup(index_row(j),2) = y(round(end+1 + f_measu_coup_row(j)/f_tick)) - Comp_O_ax2(index_row(j),1);
                            end
                        end
                    
                    end
                    
            end
        end

    otherwise
        disp("扫频模式违法");
        return
end
    


%% 处理结果，取名并存储

% 先取名。
switch INJ_MODE
    case ModeTest.P1
        name_inj = "P1";
    case ModeTest.P3_ab
        name_inj = "P3ab";
    case ModeTest.P3_dq
        name_inj = "P3dq";
    case ModeTest.P3_pn
        name_inj = "P3pn";
end

switch MES_MODE
    case ModeTest.P1
        name_mes = "P1";
    case ModeTest.P3_ab
        name_mes = "P3ab";
    case ModeTest.P3_dq
        name_mes = "P3dq";
    case ModeTest.P3_pn
        name_mes = "P3pn";
end

if isempty(save_name)
    tempt = datetime;
    save_name=string("SweepResult_"+...
        yyyymmdd(tempt)+"_"+...
        hour(tempt)+"_"+minute(tempt)+"_"+...
        round(second(tempt)));
end
save_name_last = save_name + "_" + name_inj + "_to_" + name_mes + ".mat"; % 双引号表str，单引号表char


% 由于 F_FREQ 是按照 f_freqs_pos; f_freqs_neg 的顺序排的，所以可以先拆分所测得小信号
Comp_i_pos = Comp_i(1:length(f_freqs_pos),:);
Comp_o_self_pos = Comp_o_self(1:length(f_freqs_pos),:);
Comp_o_coup_pos = Comp_o_coup(1:length(f_freqs_pos),:);
if ~isempty(f_freqs_neg)
    Comp_i_neg = Comp_i(length(f_freqs_pos)+1:end,:);
    Comp_o_self_neg = Comp_o_self(length(f_freqs_pos)+1:end,:);
    Comp_o_coup_neg = Comp_o_coup(length(f_freqs_pos)+1:end,:);
else
    Comp_i_neg = [];
    Comp_o_self_neg = [];
    Comp_o_coup_neg = [];
end

%{
% 给程序维护人员看的：截止到这里，我们有了这些存储：
% Comp_i_pos, Comp_i_neg, 如果注入是 P3_ab 或 P3_dq, 会有第2列；如果注入是 P3_pn, 会有neg
% Comp_o_self_pos, Comp_o_self_neg, 解释同上
% Comp_o_coup_pos, Comp_o_coup_neg, 解释同上；另外，仅当测量是 P3_ab 或 P3_dq 时, 它存在
% 之所以要分 pos/neg 而不是同样放到第二列，是因为 ab/dq 仅仅是同频不同轴，但 pn 正负半轴对应的是完全不同的频率，但如果放到同一个轴里又画不了bode图（bode图横坐标不能跨零）
% 但是在存储时，由于 pn/dq 建模是并列的建模形式，所以还是统一存成cell的形式，cell里分别是各响应的 frequency, magnitude, phase, 每个占一列
%}

if INJ_MODE == ModeTest.P1 && MES_MODE == ModeTest.P1
    % 注入 P1, 测量 P1
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    
    sweep_result = ...
    {"self_self"; ...
    [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180]};    % 保证相位结果在(-180,+180)之内

elseif INJ_MODE == ModeTest.P1 && (MES_MODE == ModeTest.P3_ab || MES_MODE == ModeTest.P3_dq || MES_MODE == ModeTest.P3_pn)
    % 注入 P1, 测量 P3_ab/dq 或 P3_pn
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_SC = Comp_o_coup_pos(:,1)./Comp_i_pos(:,1);
    
    sweep_result = ...
    {"self_self","self_coup"; ...
     [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_SC), mod(angle(Comp_T_SC)*180/pi + 180, 360)-180]};
    
elseif (INJ_MODE == ModeTest.P3_ab || INJ_MODE == ModeTest.P3_dq) && MES_MODE == ModeTest.P1
    % 注入 P3_ab/dq, 测量 P1
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_CS = Comp_o_self_pos(:,2)./Comp_i_pos(:,2);
    
    sweep_result = ...
    {"self_self","coup_self"; ...
     [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_CS), mod(angle(Comp_T_CS)*180/pi + 180, 360)-180]};
     
elseif (INJ_MODE == ModeTest.P3_ab || INJ_MODE == ModeTest.P3_dq) && (MES_MODE == ModeTest.P3_ab || MES_MODE == ModeTest.P3_dq || MES_MODE == ModeTest.P3_pn)
    % 注入 P3_ab/dq, 测量 P3_ab/dq 或 P3_pn
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_SC = Comp_o_coup_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_CS = Comp_o_self_pos(:,2)./Comp_i_pos(:,2);
    Comp_T_CC = Comp_o_coup_pos(:,2)./Comp_i_pos(:,2);
    
    sweep_result = ...
    {"self_self","self_coup","coup_self","coup_coup"; ...
     [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_SC), mod(angle(Comp_T_SC)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_CS), mod(angle(Comp_T_CS)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_CC), mod(angle(Comp_T_CC)*180/pi + 180, 360)-180]};
    
elseif INJ_MODE == ModeTest.P3_pn && MES_MODE == ModeTest.P1
    % 注入 P3_pn, 测量 P1
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_CS = Comp_o_self_neg(:,1)./Comp_i_neg(:,1);
    
    sweep_result = ...
    {"self_self","coup_self"; ...
     [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180], ...
     [f_freqs_neg, abs(Comp_T_CS), mod(angle(Comp_T_CS)*180/pi + 180, 360)-180]};
    
elseif INJ_MODE == ModeTest.P3_pn && (MES_MODE == ModeTest.P3_ab || MES_MODE == ModeTest.P3_dq)
    % 注入 P3_pn, 测量 P3_ab/dq
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_SC = Comp_o_coup_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_CS = Comp_o_self_neg(:,1)./Comp_i_neg(:,1);
    Comp_T_CC = Comp_o_coup_neg(:,1)./Comp_i_neg(:,1);
    
    sweep_result = ...
    {"self_self","self_coup","coup_self","coup_coup"; ...
     [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_SC), mod(angle(Comp_T_SC)*180/pi + 180, 360)-180], ...
     [f_freqs_neg, abs(Comp_T_CS), mod(angle(Comp_T_CS)*180/pi + 180, 360)-180], ...
     [f_freqs_neg, abs(Comp_T_CC), mod(angle(Comp_T_CC)*180/pi + 180, 360)-180]};
        
elseif INJ_MODE == ModeTest.P3_pn && MES_MODE == ModeTest.P3_pn
    % 注入 P3_pn, 测量 P3_pn
    Comp_T_SS = Comp_o_self_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_SC = Comp_o_coup_pos(:,1)./Comp_i_pos(:,1);
    Comp_T_CS = Comp_o_coup_neg(:,1)./Comp_i_neg(:,1);
    Comp_T_CC = Comp_o_self_neg(:,1)./Comp_i_neg(:,1);
    
    sweep_result = ...
    {"self_self","self_coup","coup_self","coup_coup"; ...
     [f_freqs_pos, abs(Comp_T_SS), mod(angle(Comp_T_SS)*180/pi + 180, 360)-180], ...
     [f_freqs_pos, abs(Comp_T_SC), mod(angle(Comp_T_SC)*180/pi + 180, 360)-180], ...
     [f_freqs_neg, abs(Comp_T_CS), mod(angle(Comp_T_CS)*180/pi + 180, 360)-180], ...
     [f_freqs_neg, abs(Comp_T_CC), mod(angle(Comp_T_CC)*180/pi + 180, 360)-180]};
    
else
    disp("没有设计这种响应测量的换算情况。请联系程序员。")
    
end
% 注意幅值是以绝对值的形式给的，不是dB

save(save_name_last, "sweep_result");