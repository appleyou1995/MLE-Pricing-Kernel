function R_grid = generate_non_uniform_grid(R_min, R_max, Target_Points)

% GENERATE_NON_UNIFORM_GRID 產生中間加密的非均勻網格
%
% 用途：
%   針對資產定價模型，在 At-the-Money (R=1) 附近加密，
%   以捕捉 Pricing Kernel 劇烈變化的部分，並精確抓出單調性違規點。
%
% 輸入：
%   R_min, R_max  : 網格的全域範圍 (例如 0.3 ~ 3.0)
%   Target_Points : 期望的總點數 (例如 10000)
%
% 輸出：
%   R_grid        : (N x 1) 的向量，已排序且去重


    % --- 設定加密區間與密度參數 ---
    % 根據經驗，Singularity 常發生在 0.8 ~ 1.2 之間
    R_focus_low  = 0.8; 
    R_focus_high = 1.2; 
    
    % 設定分配比例：讓 60% 的點都集中在核心區間
    ratio_focus = 0.60; 

    % --- 計算各區段點數 ---
    N_focus = round(Target_Points * ratio_focus);
    N_tails = Target_Points - N_focus;
    N_left  = floor(N_tails / 2);  % 左尾點數
    N_right = N_tails - N_left;    % 右尾點數

    % --- 產生三段網格 ---
    % 1. 左尾 (疏)
    if R_min < R_focus_low
        R_part1 = linspace(R_min, R_focus_low, N_left);
    else
        R_part1 = [];
    end

    % 2. 核心 (密)
    R_part2 = linspace(R_focus_low, R_focus_high, N_focus);

    % 3. 右尾 (疏)
    if R_max > R_focus_high
        R_part3 = linspace(R_focus_high, R_max, N_right);
    else
        R_part3 = [];
    end

    % --- 合併、去重、轉置 ---
    % unique 會自動排序並去除接縫處重複的點
    R_grid = unique([R_part1, R_part2, R_part3])'; 
    
end