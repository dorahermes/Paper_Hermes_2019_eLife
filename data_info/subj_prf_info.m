function [v_area,roi_labels] = subj_prf_info(subj_nr,elec)
%
% Function loads fits from pRF experiment from Jon Winawer and visual
% cortex area assigned from anatomy
%
% Output: x, y and size in degrees of visual angle 
%
% Example:
% subj_nr = 19;
% elec = 109;
% [v_area,roi_labels] = get_prf_info(subj_nr,elec);   
%
% DH 2018

roi_labels = {'?','V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2'};

switch subj_nr
        
    case 19
        switch elec
            case 106
                v_area = 6;% V2 or hV4 - ventral - foveal confluence...
            case 107
                v_area = 4;%v3 ventral
            case 108
                v_area = [2 3];%v1/v2 ventral
            case 109
                v_area = 2;%v1 ventral
            case 115
                v_area = [2 3];%v1/v2 dorsal
            case 119
                v_area = 3;%v2 dorsal
            case 120
                v_area = [3 4];%v2/v3 dorsal
            case 121
                v_area = [4];%v3 dorsal
            case 122
                v_area = [5];%V3B dorsal
        end
        
    case 24
        switch elec
            case 45
                v_area = 2;%v1 ventral
            case 46
                v_area = 2;%v1 ventral
        end
        
     case 1001
        switch elec
            case 37
                v_area = 5;%v3a dorsal, at the upper edge
            case 49
                v_area = 3;%v2 ventral
            case 50
                v_area = 3;%v2 ventral
            case 51
                v_area = 2;%v1 dorsal
            case 52
                v_area = [2 3];%v1/2 ventral
            case 57
                v_area = [3];%v2 ventral
            case 58
                v_area = [2];%v1 ventral
            case 59
                v_area = [2];%v1 ventral
            case 60
                v_area = [2 3];%v1/2 ventral
        end
end

    