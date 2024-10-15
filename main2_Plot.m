%% simulink 扫频程序后续：画图
% 
% 因为画图需求往往与扫频不一致，所以把扫频主程序截止到存储为止，用这个程序读取与画图。

%% 设置部分

%{
% 如果你有传函，要和扫频结果一起显示，可以把它存成类似
%   sys_cell = {"self_self", "self_coup"; tf('s'), tf([0,1],[1,0])};
%   save("my_sys.mat","sys_cell");
% 的形式。然后把下面这个设置成 true ，来提醒程序读取。
%}
TF_PLOT = false;
TF_name = "my_sys";
TF_plot_option = "k-";  % 统一设置的画图选项
swp_plot_option = "r*";

load_name = "Swp_A1_P1_to_P1";   % 读取之前存储的扫频结果

% 测量响应的时候，我存储的是绝对值。如果你希望用dB来画，那么请把这句设为true。
PLOT_IN_DB = true;

%% 画图部分

% 准备画布
FontSize = 8;

gcf_bode = figure(2); clf; set(gcf,"Units","centimeters","Position",[2,6,16,12])
gca_handles = cell(2,4);
line_handles = cell(2,4);

% 最先读扫频结果。就算要画传函，也需要扫频结果提供的频点集。
load(load_name + ".mat")
[m_res, n_res] = size(sweep_result);

% 先画传函（如果有）
if TF_PLOT
    load(TF_name + ".mat")
    [m_sys, n_sys] = size(sys_cell);

    for i = 1:n_sys
        % 共通的定义子图部分
        switch sweep_result{1, i}
            case "self_self"
                if ~(exist("gca_mag_SS","var") && isvalid(gca_mag_SS))
                    gca_mag_SS = subplot("Position", [0.06, 0.75, 0.42, 0.22]);
                    gca_pha_SS = subplot("Position", [0.06, 0.53, 0.42, 0.22]);
                    gca_handles{1,1} = gca_mag_SS;
                    gca_handles{2,1} = gca_pha_SS;
                end
                gca_mag = gca_mag_SS;
                gca_pha = gca_pha_SS;
            case "self_coup"
                if ~(exist("gca_mag_SC","var") && isvalid(gca_mag_SC))
                    gca_mag_SC = subplot("Position", [0.56, 0.75, 0.42, 0.22]);
                    gca_pha_SC = subplot("Position", [0.56, 0.53, 0.42, 0.22]);
                    gca_handles{1,2} = gca_mag_SC;
                    gca_handles{2,2} = gca_pha_SC;
                end
                gca_mag = gca_mag_SC;
                gca_pha = gca_pha_SC;
            case "coup_self"
                if ~(exist("gca_mag_CS","var") && isvalid(gca_mag_cS))
                    gca_mag_CS = subplot("Position", [0.06, 0.27, 0.42, 0.22]);
                    gca_pha_CS = subplot("Position", [0.06, 0.05, 0.42, 0.22]);
                    gca_handles{1,3} = gca_mag_CS;
                    gca_handles{2,3} = gca_pha_CS;
                end
                gca_mag = gca_mag_CS;
                gca_pha = gca_pha_CS;
            case "coup_coup"
                if ~(exist("gca_mag_CC","var") && isvalid(gca_mag_CC))
                    gca_mag_CC = subplot("Position", [0.56, 0.27, 0.42, 0.22]);
                    gca_pha_CC = subplot("Position", [0.56, 0.05, 0.42, 0.22]);
                    gca_handles{1,4} = gca_mag_CC;
                    gca_handles{2,4} = gca_pha_CC;
                end
                gca_mag = gca_mag_CC;
                gca_pha = gca_pha_CC;
            otherwise
                disp("第"+ i +"组传函设置无效")
                return;
        end

        set(gca_mag, "NextPlot", "add");
        set(gca_pha, "NextPlot", "add");
        
        % 计算传函等效的响应集
        freq_sys = [];
        mag_sys = [];
        pha_sys = [];
        for j = 1:n_res
            if strcmp(sys_cell{1, i}, sweep_result{1, j})
                freq_sys = sweep_result{2, j}(:,1);
                temp_res = arrayfun(@(w) evalfr(sys_cell{2, i}, 1j*w), freq_sys*2*pi);
                mag_sys = abs(temp_res);
                pha_sys = mod(angle(temp_res)*180/pi + 180, 360)-180;
                break
            end
        end
        clear temp_res
        
        freq_abs = abs(freq_sys);

        % 以下是共通的画图命令
        if PLOT_IN_DB
            line_mag = semilogx(gca_mag,freq_abs,20*log10(mag_sys),TF_plot_option);
        else
            line_mag = semilogx(gca_mag,freq_abs,mag_sys,TF_plot_option);
        end
        line_pha = semilogx(gca_pha,freq_abs,pha_sys,TF_plot_option);
        line_handles{1,end+1} = line_mag;
        line_handles{2,end+1} = line_pha;

        if PLOT_IN_DB
            myBodeFix(gca_mag,gca_pha,min(freq_abs),max(freq_abs),FontSize,"dB")
        else
            myBodeFix(gca_mag,gca_pha,min(freq_abs),max(freq_abs),FontSize,"abs")
        end

        if freq_sys(1) < 0
            containsDash = cellfun(@(x) any(x == '-'), gca_pha.XTickLabel);
            if ~any(containsDash)
                set(gca_pha, "XTickLabel", cellfun(@(x) ['-', x], gca_pha.XTickLabel, 'UniformOutput', false));
            end
            set(gca_mag, "XDir", "reverse");
            set(gca_pha, "XDir", "reverse");
        end
    end
end

% 再画扫频点集
for j =1:n_res
    % 共通的定义子图部分
    switch sweep_result{1, j}
        case "self_self"
            if ~(exist("gca_mag_SS","var") && isvalid(gca_mag_SS))
                gca_mag_SS = subplot("Position", [0.06, 0.75, 0.42, 0.22]);
                gca_pha_SS = subplot("Position", [0.06, 0.53, 0.42, 0.22]);
                gca_handles{1,1} = gca_mag_SS;
                gca_handles{2,1} = gca_pha_SS;
            end
            gca_mag = gca_mag_SS;
            gca_pha = gca_pha_SS;
        case "self_coup"
            if ~(exist("gca_mag_SC","var") && isvalid(gca_mag_SC))
                gca_mag_SC = subplot("Position", [0.56, 0.75, 0.42, 0.22]);
                gca_pha_SC = subplot("Position", [0.56, 0.53, 0.42, 0.22]);
                gca_handles{1,2} = gca_mag_SC;
                gca_handles{2,2} = gca_pha_SC;
            end
            gca_mag = gca_mag_SC;
            gca_pha = gca_pha_SC;
        case "coup_self"
            if ~(exist("gca_mag_CS","var") && isvalid(gca_mag_CS))
                gca_mag_CS = subplot("Position", [0.06, 0.27, 0.42, 0.22]);
                gca_pha_CS = subplot("Position", [0.06, 0.05, 0.42, 0.22]);
                gca_handles{1,3} = gca_mag_CS;
                gca_handles{2,3} = gca_pha_CS;
            end
            gca_mag = gca_mag_CS;
            gca_pha = gca_pha_CS;
        case "coup_coup"
            if ~(exist("gca_mag_CC", "var") && isvalid(gca_mag_CC))
                gca_mag_CC = subplot("Position", [0.56, 0.27, 0.42, 0.22]);
                gca_pha_CC = subplot("Position", [0.56, 0.05, 0.42, 0.22]);
                gca_handles{1,4} = gca_mag_CC;
                gca_handles{2,4} = gca_pha_CC;
            end
            gca_mag = gca_mag_CC;
            gca_pha = gca_pha_CC;
        otherwise
            disp("第"+ j +"组数据设置无效")
            return;
    end
    
    set(gca_mag, "NextPlot", "add");
    set(gca_pha, "NextPlot", "add");
    
    % 读取响应集
    freq_T = sweep_result{2, j}(:,1); mag_T = sweep_result{2, j}(:,2); pha_T = sweep_result{2, j}(:,3);
    freq_abs = abs(freq_T);  % 横轴频率先不管正负，按绝对值画图

    % 以下是共通的画图命令
    if PLOT_IN_DB
        line_mag = semilogx(gca_mag,freq_abs,20*log10(mag_T),"r*");
    else
        line_mag = semilogx(gca_mag,freq_abs,mag_T,"r*");
    end
    line_pha = semilogx(gca_pha,freq_abs,pha_T,"r*");
    line_handles{1,end+1} = line_mag;
    line_handles{2,end+1} = line_pha;
    
    if PLOT_IN_DB
        myBodeFix(gca_mag,gca_pha,min(freq_abs),max(freq_abs),FontSize,"dB")
    else
        myBodeFix(gca_mag,gca_pha,min(freq_abs),max(freq_abs),FontSize,"abs")
    end

    if freq_T(1) < 0    % 补充对横轴频率为负时的补充标注与坐标轴反转
        containsDash = cellfun(@(x) any(x == '-'), gca_pha.XTickLabel);
        if ~any(containsDash)
            set(gca_pha, "XTickLabel", cellfun(@(x) ['-', x], gca_pha.XTickLabel, 'UniformOutput', false));
        end
        set(gca_mag, "XDir", "reverse");
        set(gca_pha, "XDir", "reverse");
    end
end

