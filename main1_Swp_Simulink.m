%% simulink 扫频程序 支持 SIMULINK 与 PLECS
% 
% 使用说明（减缩版）：
%   plecs模型用4.7以上的版本打开挂在后台。建议用变步长仿真。
%   在 PORT_MODE = "AC" 下，输出有7个维度： inA, inB, inC, outA, outB, outC, theta. 不用的信号请用0占位。
%   在 PORT_MODE = "DC" 下，输出有2个维度： in, out. 这个模式是为了减小数据传输量做的。
%   请记得把模型总运行时间设为runT，并且不要 在 initialization 里重新定义以 inj 开头的变量。不过，runT 可以重定义。
%   仅支持单相扫频，且不支持一次性扫多个频点。如果你想扫dq，那请在仿真模型的输入输出中做变换。

%% 扫频设置
% 如果你是用户，请你只改动这部分的变量，没有特殊需求就不要改其它部分了。
sim_mode = "PLECS";   % PLECS or SIMULINK
path = pwd;     % 在当前路径定位plecs模型
model_name = 'VENA_sim_2';    % 模型路径

signal_introduction = "电压闭环，前馈带2滤波";    % 本次测量的响应的备注
plot_option = "-b*";    % 本次测量的响应的画图标注
save_name = "output_sim_swp/AFE_Vin2Iin_2";    % 扫频结果保存名。若缺省为""，则按当前时间保存。后面保存的时候会自动加上模式注释。

% 如果还需要调用别的路径下的自定义函数，请在这里 addpath
addpath(genpath('.\Includes'));   % 添加函数路径 

% 下述“注入端口”不必须是仿真中真正注入谐波的端口。但在三相双轴测试模式下，需要保证在谐波源单轴注入时，“注入端口”的谐波没有混杂。
INJ_MODE = ModeTest.P1;     % 注入端口测量模式
MES_MODE = ModeTest.P1;     % 测量端口测量模式
PORT_MODE = "AC";  % "AC" or "DC"

COUP_MODE = 0;   % 测量端口相对注入端口的基准频率耦合情况。
% 即，测量多端口频率耦合时，测量端口使用的 F_FREQS' 是否需要做特殊位移。 
% -2, -1, 0, 1, 2 -- "0"表示同样频率。在此基础上，每个数表示测量端口额外加减的基波频率个数。即， F_FREQS' = F_FREQS + COUP_MODE * F_FUND。
% （备用功能）98, 99, 100, 101, 102 -- "100"表示测量端口是注入端口关于基波的对称。在此基础上，加减基波个数。即， F_FREQS' = (2*F_FUND-F_FREQS) + (COUP_MODE-100) * F_FUND。

AUTO_MODE = false;  % AUTO模式。若启动该选项，则认为注入、测量的端口性质相同，自动设置 MES_MODE = INJ_MODE, COUP_MODE = 0。

injA_ref = 10;   % 注入小信号幅值 ###### 重要 ######

% round( logspace( log10(最小频率),log10(最大频率),取点数量 )/频率分辨率 )*频率分辨率
F_FREQS = unique([...
    round( logspace(log10(1),   log10(40), 14)/1  )*1, ...
    round( logspace(log10(40),   log10(1e2), 6)/5  )*5, ...
    round( logspace(log10(1e2), log10(1e3), 10)/25  )*25, ...
    round( logspace(log10(1e3), log10(1e4), 10)/50  )*50])';
% F_FREQS = unique([...
%     round( logspace(log10(1e1), log10(1e2), 10)/10  )*10, ...
%     round( logspace(log10(1e2), log10(1e3), 10)/20  )*20, ...
%     round( logspace(log10(1e3), log10(1e4), 10)/50  )*50])';

% 该频率是注入端口的频率。注意请把它设置为竖的向量。
% 在dq模式下，它是在dq系下的等效注入频率。
% 如果你想注入负序，请直接用负频率代替。若它同时包含正频率与负频率，下面画图时会分出一个正半图，一个负半图。
% 在正负序，AUTO模式下，它会取该频率的绝对值作为正序注入集，并自动生成负序。
% figure(1); semilogx(F_FREQS,ones(1,length(F_FREQS)),"*")   % 画个图验证扫频点排布，可以注释掉

F_TICK = double(gcd(sym(F_FREQS)));    % 频率分辨率。它防止高频测试中忽略基波和纹波的影响。严格来说，请把它设成基波频率的因数。
F_FUND = 50;    % 基波频率
T_STEADY = 0.1;   % 认为从此开始，仿真模型处于稳态
T_SIM_MIN = 1e-7;   % 仿真中的时间最小值（该扫频程序允许离散与连续仿真）

%% 对仿真设置的处理

% 去零
F_FREQS(F_FREQS == 0) = [];
f_freqs_pos = sort(F_FREQS(F_FREQS>0));
f_freqs_neg = -sort(-(F_FREQS(F_FREQS < 0)));
F_FREQS = [f_freqs_pos;f_freqs_neg];    % 排版：先是从小到大排序的正数部分，再是绝对值从小到大排序的负数部分

N_FREQS = length(F_FREQS);

% 由于 F_FREQS 可能发生了变化，所以重新排版（排版信息在后续换算结果时会使用到）
f_freqs_pos = sort(F_FREQS(F_FREQS>0));
f_freqs_neg = -sort(-(F_FREQS(F_FREQS < 0)));
F_FREQS = [f_freqs_pos;f_freqs_neg];

% 生成实际的测量频点（中心频点）
if COUP_MODE < 50
    F_MEASU = F_FREQS + COUP_MODE * F_FUND;
else
    F_MEASU = (2*F_FUND-F_FREQS) + (COUP_MODE-100) * F_FUND;
end

F_MEASU_COUP = 2*F_FUND-F_MEASU;    % 仅为 P3_pn 模式的测量准备

% 去掉无效频点（F_FREQS, F_MEASU, F_MEASU_COUP 中任何一个可能为0的项）
F_MEASU_COUP(F_FREQS == 0) = [];
F_MEASU(F_FREQS == 0) = [];
F_FREQS(F_FREQS == 0) = [];
F_MEASU_COUP(F_MEASU == 0) = [];
F_FREQS(F_MEASU == 0) = [];
F_MEASU(F_MEASU == 0) = [];
if ~(numel(F_FREQS) == numel(F_MEASU) && numel(F_FREQS) == numel(F_MEASU_COUP))
    disp('频率去重部分出错，请联系调试')
    return
end

k_inj = F_FREQS/F_FUND;
k_mes = F_MEASU/F_FUND;

% 生成最大的测量频点
f_mes_abs_max = max(max(F_FREQS), -min(F_FREQS));
f_temp = [F_MEASU; F_MEASU_COUP];
f_mes_abs_max = max([f_mes_abs_max, max(f_temp), -min(f_temp)]);

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

% 占缓存
Comp_i = zeros(N_FREQS, 1);
Comp_I = zeros(N_FREQS, 1);
Comp_o_self = zeros(N_FREQS, 1);
Comp_o_coup = zeros(N_FREQS, 1);


Comp_O_ax1 = zeros(N_FREQS, 1);
Comp_O_ax2 = zeros(N_FREQS, 1);

fprintf('\n共需测试%d次\n', N_FREQS);

% 建线程，开模型
if(sim_mode == "PLECS")
    proxy = jsonrpc('http://localhost:1080', 'Timeout', 300);
    proxy.plecs.load([path '/' model_name '.plecs']);
    % proxy.plecs.scope([model_name '/Scope'], 'ClearTraces');
end

%% 第一次运行仿真：获取稳态工作点矢量
T_SCALE = 1/F_TICK;     % 根据取频点的分辨率计算需要取的时间段长度
% T_SCALE = round(1/K_PER)/F_FUND;

% Signal 模式，构造设置并仿真
if(sim_mode == "SIMULINK")
    simIn = Simulink.SimulationInput(model_name);
    simIn = setVariable(simIn, "runT", T_STEADY+T_SCALE);
    simIn = setVariable(simIn, "injA", 0);
    simIn = setVariable(simIn, "injF", 50);
    simIn = setVariable(simIn, "injM", int32(ModeTest.signal));
    simOut = sim(simIn);
    result_time = simOut.freq_swp.time; % 列向量
    result_signal = simOut.freq_swp.signals.values; % 每列是一个变量，共7列
elseif(sim_mode == "PLECS")
    simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', 50, 'injA', 0));
    result = proxy.plecs.simulate(model_name, simStruct);   % 运行仿真
    result_time = transpose(result.Time);
    result_signal = transpose(result.Values); % 每行是一个变量，共7行。转置一下统一格式
end
[result_m, result_n] = size(result_signal);      % 用于验证维度是否与模式设置匹配
% if result_n < 7
%     msgbox("输出通道数不足！");
%     return
% end

% 不管仿真是定步长、变步长，对采样线性插值就行
t = T_STEADY : T_TICK : result_time(end);      N_t = length(t);    % 取稳态时间段

switch INJ_MODE
    case ModeTest.P1
        % 获取时域测量，采样插值以模拟定点示波器
        x = interp1(result_time, result_signal(:,1), t);
        y = fft(x)*2/N_t;
        
        % 提取响应
        for i=1:N_FREQS
            if F_FREQS(i) >= 0
                Comp_I(i) = y(round(1 + F_FREQS(i)/F_TICK));
            else
                Comp_I(i) = y(round(end+1 + F_FREQS(i)/F_TICK));
            end
        end
end

switch MES_MODE
    case ModeTest.P1
        % 获取时域测量，采样插值以模拟定点示波器
        if PORT_MODE == "DC"
            x = interp1(result_time, result_signal(:,2), t);
        else
            x = interp1(result_time, result_signal(:,4), t);
        end
        y = fft(x)*2/N_t;
        
        % 提取响应
        for i=1:length(F_MEASU)
            if F_MEASU(i) >= 0
                Comp_O_ax1(i) = y(round(1 + F_MEASU(i)/F_TICK));
            else
                Comp_O_ax1(i) = y(round(end+1 + F_MEASU(i)/F_TICK));
            end
        end
end

fprintf('已获得稳态基准\n');

%% 遍历运行仿真：测量小信号响应

for i = 1:N_FREQS
    f_freq = F_FREQS(i, :);
    f_meas = F_MEASU(i);
    f_meas_coup = F_MEASU_COUP(i);

    % 计算当前注入行对应的分辨率
    f_tick = F_TICK;
    if(mod(f_freq, F_FUND)==0 && mod(f_meas, F_FUND)==0 && mod(f_meas_coup, F_FUND)==0)
        f_tick = F_FUND;
    elseif(mod(f_freq, 10)==0 && mod(f_meas, 10)==0 && mod(f_meas_coup, 10)==0)
        f_tick = 10; 
    end
    T_SCALE = 1/f_tick;

    % 决定当前注入行的模式
    switch INJ_MODE
        case ModeTest.P1
            injM_ref = ModeTest.P1;
    end

    % 构造设置并仿真
    if(sim_mode == "SIMULINK")
        simIn = Simulink.SimulationInput(model_name);
        simIn = setVariable(simIn, "runT", T_STEADY+T_SCALE);
        simIn = setVariable(simIn, "injA", injA_ref);
        simIn = setVariable(simIn, "injF", f_freq);
        simIn = setVariable(simIn, "injM", int32(injM_ref));
        simOut = sim(simIn);
        result_time = simOut.freq_swp.time; % 列向量
        result_signal = simOut.freq_swp.signals.values; % 每列是一个变量，共7列
    elseif(sim_mode == "PLECS")
        simStruct = struct('ModelVars', struct('runT', T_STEADY+T_SCALE, 'injF', f_freq, 'injA', injA_ref));
        result = proxy.plecs.simulate(model_name, simStruct);   % 运行仿真
        result_time = transpose(result.Time);
        result_signal = transpose(result.Values); % 每行是一个变量，共7行。转置一下统一格式
    end

    % 提取结果
    t = T_STEADY : T_TICK : result_time(end);      N_t = length(t);

    switch INJ_MODE
        case ModeTest.P1

            x = interp1(result_time, result_signal(:,1), t);
            y = fft(x)*2/N_t;
            
            % 提取响应
            if f_freq > 0
                Comp_i(i) = y(round(1 + f_freq/f_tick)) - Comp_I(i);
            else
                Comp_i(i) = y(round(end+1 + f_freq/f_tick)) - Comp_I(i);
            end

    end

    switch MES_MODE
        case ModeTest.P1

            % 获取时域测量，采样插值以模拟定点示波器
            if PORT_MODE == "DC"
                x = interp1(result_time, result_signal(:,2), t);
            else
                x = interp1(result_time, result_signal(:,4), t);
            end
            y = fft(x)*2/N_t;
            
            % 提取响应

            if f_meas > 0
                Comp_o_self(i) = y(round(1 + f_meas/f_tick)) - Comp_O_ax1(i);
            else
                Comp_o_self(i) = y(round(end+1 + f_meas/f_tick)) - Comp_O_ax1(i);
            end

    end

    temp_resp = abs(Comp_o_self(i)) / abs(Comp_i(i));
    fprintf('第%d次测试已结束，测得幅值为 %.2f dB\n', [i, 20*log( temp_resp )/log(10)]);
end

    


%% 处理结果，取名并存储

% 先取名。
if isempty(save_name) || strcmp(save_name,"")
    tempt = datetime;
    save_name=string("SwpRes_"+...
        yyyymmdd(tempt)+"_"+...
        hour(tempt)+"_"+minute(tempt)+"_"+...
        round(second(tempt)));
end

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

end
% 注意幅值是以绝对值的形式给的，不是dB

sweep_result{3,1} = plot_option;
sweep_result{4,1} = signal_introduction;
save(save_name + ".mat", "sweep_result");

msgbox("你扫频跑完了！");
