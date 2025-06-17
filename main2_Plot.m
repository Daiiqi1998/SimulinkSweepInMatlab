%% simulink 扫频程序后续：画图
% 
% 因为画图需求往往与扫频不一致，所以把扫频主程序截止到存储为止，用这个程序读取与画图。

%% 设置部分

%{
% 如果你有传函，可以把它存为类似
%   tf_result = {"目标阻抗"; tf('s'); "b-"};
%   save("my_sys.mat","tf_result");
% 的形式，然后在运行该画图程序前双击读取。
% 当有多个传函时，每个传函占一列。
%}
TF_PLOT = true;

%% 画图部分

% 准备画布
FontSize = 8;

gcf_bode = figure(2); clf; set(gcf,"Units","centimeters","Position",[2,6,16,12])
gca_handles = cell(2,4);
line_handles = cell(2,4);

% 最先读扫频结果。就算要画传函，也需要扫频结果提供的频点集。
[m_swp, n_swp] = size(sweep_result);

%% 再画扫频点集
for j =1:n_swp
    % 共通的定义子图部分
    if ~(exist("gca_mag","var") && isvalid(gca_mag))
        gca_mag = subplot("Position", [0.06, 0.75, 0.42, 0.22]);
        gca_pha = subplot("Position", [0.06, 0.53, 0.42, 0.22]);
    end
    set(gca_mag, "NextPlot", "add");
    set(gca_pha, "NextPlot", "add");
    
    % 读取响应集
    freq_T = sweep_result{2, j}(:,1); mag_T = sweep_result{2, j}(:,2); pha_T = sweep_result{2, j}(:,3);
    freq_abs = abs(freq_T);  % 横轴频率先不管正负，按绝对值画图

    % 以下是共通的画图命令
    line_mag = semilogx(gca_mag,freq_abs,20*log10(mag_T),sweep_result{3, j});

    line_pha = semilogx(gca_pha,freq_abs,pha_T,sweep_result{3, j});
    line_handles{1,end+1} = line_mag;
    line_handles{2,end+1} = line_pha;
    
    myBodeFix(gca_mag,gca_pha,min(freq_abs),max(freq_abs),FontSize,"dB")
	
	% 相位图lim设置
    set(gca_pha,'YLim',[-180 180]);
	set(gca_pha,'YTick',[-180 -90 0 90 180]);
end

%% 先画传函（如果有）
if TF_PLOT
    [m_tf, n_tf] = size(tf_result);

    for i = 1:n_tf
        % 共通的定义子图部分
        if ~(exist("gca_mag","var") && isvalid(gca_mag))
            gca_mag = subplot("Position", [0.06, 0.75, 0.42, 0.22]);
            gca_pha = subplot("Position", [0.06, 0.43, 0.42, 0.22]);
        end
        set(gca_mag, "NextPlot", "add");
        set(gca_pha, "NextPlot", "add");
        
        % 计算传函等效的响应集
        % freq_sys = sweep_result{2, 1}(:,1); freq_sys = [freq_sys; (freq_sys(1:end-1) + freq_sys(2:end)) / 2];  freq_sys = unique(freq_sys(:));
        freq_sys = logspace(1, 4, 400); % logspace 里面填 10 的指数
        
        [mag_sys, pha_sys] = bode(tf_result{2, i}, freq_sys*2*pi); mag_sys = mag_sys(1, :); pha_sys = mod(pha_sys(1, :)+180,360)-180;
        
        freq_abs = abs(freq_sys);

        % 以下是共通的画图命令

        line_mag = semilogx(gca_mag,freq_abs,20*log10(mag_sys),tf_result{3, i});

        line_pha = semilogx(gca_pha,freq_abs,pha_sys,tf_result{3, i});
        line_handles{1,end+1} = line_mag;
        line_handles{2,end+1} = line_pha;

	    % 相位图lim设置
        % set(gca_pha,'YLim',[-180 180]);
	    % set(gca_pha,'YTick',[-180 -90 0 90 180]);

        % set(gca_pha,'XLim',[1 1e5]);
        % set(gca_mag,'XLim',[1 1e5]);
        % set(gca_pha,'XTick',[1 10 1e2 1e3 1e4 1e5]);
        % set(gca_mag,'XTick',[1 10 1e2 1e3 1e4 1e5]);
        
        set(gca_pha,'YLim',[-180 90]);
	    set(gca_pha,'YTick',[-180 -90 0 90]);
    end
end

