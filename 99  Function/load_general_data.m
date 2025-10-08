function [Risk_Free_Rate, Smooth_AllR, Smooth_AllR_RND] = ...
    load_general_data(varargin)

    Path_Data = 'D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method';

    Target_TTM = 30;
    
    % Load risk-free rate R_f^t
    Path_Data_01_main = fullfile(Path_Data, 'Code', '01  原始資料處理');
    FileName = 'Risk_Free_Rate.csv';
    Risk_Free_Rate_All = readtable(fullfile(Path_Data_01_main, FileName));
    Risk_Free_Rate = Risk_Free_Rate_All.rf_gross_29d;
    
    % Load Q-measure PDF tables: R axis and corresponding f^*_t(R)
    Path_Data_02 = fullfile(Path_Data, 'Code', '02  輸出資料');
    Smooth_AllR = [];
    Smooth_AllR_RND = [];
    
    years_to_merge = 1996:2021;
    for year = years_to_merge
        input_filename = fullfile(Path_Data_02, sprintf('TTM_%d_RND_Tables_%d.mat', Target_TTM, year));
        if exist(input_filename, 'file')
            data = load(input_filename);
            Smooth_AllR = [Smooth_AllR, data.Table_Smooth_AllR];
            Smooth_AllR_RND = [Smooth_AllR_RND, data.Table_Smooth_AllR_RND];
        else
            warning('File %s does not exist.', input_filename);
        end
    end
end