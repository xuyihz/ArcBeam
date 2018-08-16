%% function
% MGT model
%
% Xu Yi, 2018

%%
function CB_model(fileID)
%% SECTION
SEC_400x400_num = '1'; SEC_300X500_num = '2';
SEC_400x400_SNAME = '400x400'; SEC_300X500_SNAME = '300X500';
SEC_400x400_H = '400'; SEC_300X500_H = '500';
SEC_400x400_B = '400'; SEC_300X500_B = '300';
SEC_TYPE = 'DBUSER'; SEC_OFFSET = 'CC'; SEC_0 = '0'; SEC_2 = '2'; SEC_YES = 'YES'; SEC_NO = 'NO'; SEC_SB = 'SB';

fprintf(fileID,'*SECTION    ; Section\n');
fprintf(fileID,'; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, [DATA1], [DATA2]                    ; 1st line - DB/USER\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, BLT, D1, ..., D8, iCEL              ; 1st line - VALUE\n;       AREA, ASy, ASz, Ixx, Iyy, Izz                                               ; 2nd line\n;       CyP, CyM, CzP, CzM, QyB, QzB, PERI_OUT, PERI_IN, Cy, Cz                     ; 3rd line\n;       Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4, Zyy, Zzz                                    ; 4th line\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, ELAST, DEN, POIS, POIC, SF, THERMAL ; 1st line - SRC\n;       D1, D2, [SRC]                                                               ; 2nd line\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, 1, DB, NAME1, NAME2, D1, D2         ; 1st line - COMBINED\n; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, 2, D11, D12, D13, D14, D15, D21, D22, D23, D24\n; iSEC, TYPE, SNAME, [OFFSET2], bSD, bWE, SHAPE, iyVAR, izVAR, STYPE                ; 1st line - TAPERED\n;       DB, NAME1, NAME2                                                            ; 2nd line(STYPE=DB)\n;       [DIM1], [DIM2]                                                              ; 2nd line(STYPE=USER)\n;       D11, D12, D13, D14, D15, D16, D17, D18                                      ; 2nd line(STYPE=VALUE)\n;       AREA1, ASy1, ASz1, Ixx1, Iyy1, Izz1                                         ; 3rd line(STYPE=VALUE)\n;       CyP1, CyM1, CzP1, CzM1, QyB1, QzB1, PERI_OUT1, PERI_IN1, Cy1, Cz1           ; 4th line(STYPE=VALUE)\n;       Y11, Y12, Y13, Y14, Z11, Z12, Z13, Z14, Zyy1, Zyy2                          ; 5th line(STYPE=VALUE)\n;       D21, D22, D23, D24, D25, D26, D27, D28                                      ; 6th line(STYPE=VALUE)\n;       AREA2, ASy2, ASz2, Ixx2, Iyy2, Izz2                                         ; 7th line(STYPE=VALUE)\n;       CyP2, CyM2, CzP2, CzM2, QyB2, QzB2, PERI_OUT2, PERI_IN2, Cy2, Cz2           ; 8th line(STYPE=VALUE)\n;       Y21, Y22, Y23, Y24, Z21, Z22, Z23, Z24, Zyy2, Zzz2                          ; 9th line(STYPE=VALUE)\n; [DATA1] : 1, DB, NAME or 2, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10\n; [DATA2] : CCSHAPE or iCEL or iN1, iN2\n; [SRC]  : 1, DB, NAME1, NAME2 or 2, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, iN1, iN2\n; [DIM1], [DIM2] : D1, D2, D3, D4, D5, D6, D7, D8\n; [OFFSET] : OFFSET, iCENT, iREF, iHORZ, HUSER, iVERT, VUSER\n; [OFFSET2]: OFFSET, iCENT, iREF, iHORZ, HUSERI, HUSERJ, iVERT, VUSERI, VUSERJ\n');

fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_400x400_num,SEC_TYPE,SEC_400x400_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_SB,SEC_2,SEC_400x400_H,SEC_400x400_B,...
    SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0);
fprintf(fileID,'   %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
    SEC_300X500_num,SEC_TYPE,SEC_300X500_SNAME,...
    SEC_OFFSET,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_YES,SEC_NO,...
    SEC_SB,SEC_2,SEC_300X500_H,SEC_300X500_B,...
    SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0,SEC_0);
fprintf(fileID,'\n');

%% append Model

%% NODE
fprintf(fileID,'*NODE    ; Nodes\n');
fprintf(fileID,'; iNO, X, Y, Z\n');

Bspan = 8400;   % 柱间跨度
levelZaxis = [0, 4000]; % 层高Z坐标
CoC = [0, 0];   % 圆心
Arc_itvl = 1000; % 定义“以直代曲”的最大直线段长度。
XYcoor = [Bspan/2,Bspan/2; -Bspan/2,Bspan/2; -Bspan/2,-Bspan/2; Bspan/2,-Bspan/2; ];    % XY坐标
lengthlevelZaxis = length(levelZaxis);  % 层数
XYcoor_num = length(XYcoor);    % 每层节点数

iNO = 0;    % 初始化
iNO_init = iNO;
for i = 1:lengthlevelZaxis	% length(A(:)) A向量元素个数
    for j = 1:XYcoor_num
        iNO = iNO+1;
        fprintf(fileID,'   %d, %.4f, %.4f, %.4f\n',...
            iNO,XYcoor(j,1),XYcoor(j,2),levelZaxis(i));
    end
end
% 以直代曲部分
iNO_main_end = iNO; % 主节点终点备份，即以直代曲节点起点备份。
XY_Deg_num = zeros(lengthlevelZaxis,XYcoor_num); % 以直代曲的各曲线的分隔节点数
P_start = zeros(1,2); P_end = zeros(1,2);
for i = 2:lengthlevelZaxis
    for j = 1:XYcoor_num %
        P_start(:) = XYcoor(j,:);
        if j == XYcoor_num
            P_end(:) = XYcoor(1,:);
        else
            P_end(:) = XYcoor(j+1,:);
        end
        [iNO, Deg_num] = CB_arc_FE(fileID, iNO, levelZaxis(i), CoC, P_start, P_end, Arc_itvl);
        XY_Deg_num(i,j) = Deg_num;
    end
end
% 以直代曲部分
iNO_end = iNO;
fprintf(fileID,'\n');

%% ELEMENT(frame) columns
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_colu = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 2; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 柱；iPRO = 1 截面编号1。
fprintf(fileID,'; 柱\n');
ELE_iPRO = 1;
iNO = iNO_init; % 初始化iNO
iEL = 0;    % 初始化
for i = 1:(lengthlevelZaxis-1)	% length(A(:)) A向量元素个数
    for j = 1:XYcoor_num
        iEL = iEL+1;
        iN1 = iNO+j;
        iN2 = iN1+XYcoor_num;
        fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
            iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
            iN1, iN2,...    % 柱单元的两个节点号
            ELE_ANGLE, ELE_iSUB);
    end
end
fprintf(fileID,'\n');

%% ELEMENT(frame) beams
fprintf(fileID,'*ELEMENT    ; Elements\n');
fprintf(fileID,'; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, iOPT(EXVAL2) ; Frame  Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element\n; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element\n');

% iEL_init_beam = iEL;
ELE_TYPE = 'BEAM'; ELE_iMAT = 2; ELE_ANGLE = 0; ELE_iSUB = 0;  % iMAT = 1材料钢结构Q345

% 横向主梁；iPRO = 2 截面编号2。
fprintf(fileID,'; 主梁\n');
ELE_iPRO = 2;
iNO = iNO_init; % 初始化iNO
iNO_arc = iNO_main_end; % 初始化
iEL_beam_0 = iEL;
for i = 1:lengthlevelZaxis	% 此行与柱单元不同，柱单元为i-1
    for j = 1:XYcoor_num
        iN1_bkp = iNO+j+XYcoor_num*(i-1);
        if j ~= XYcoor_num
            iN2_bkp = iN1_bkp+1;
        else % j = lengthXYcoor_f 时， 连接的是本环的第一个点，而不是上层内环的第一个点。
            iN2_bkp = iN1_bkp+1-XYcoor_num;
        end
        
        if XY_Deg_num(i,j) == 1 % 即此处未进行以直代曲分隔
            iEL = iEL+1;
            iN1 = iN1_bkp;
            iN2 = iN2_bkp;
            fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                iN1, iN2,...    % 梁单元的两个节点号
                ELE_ANGLE, ELE_iSUB);
        else
            for k = 1:XY_Deg_num(i,j)
                iEL = iEL+1;
                if k == 1
                    iNO_arc = iNO_arc+1;
                    iN1 = iN1_bkp;
                    iN2 = iNO_arc;
                elseif k == XY_Deg_num(i,j)
                    iN1 = iNO_arc;
                    iN2 = iN2_bkp;
                else
                    iN1 = iNO_arc;
                    iNO_arc = iNO_arc+1;
                    iN2 = iNO_arc;
                end
                fprintf(fileID,'   %d, %s, %d, %d, %d, %d, %d, %d\n',...
                    iEL, ELE_TYPE, ELE_iMAT, ELE_iPRO,...
                    iN1, iN2,...    % 梁单元的两个节点号
                    ELE_ANGLE, ELE_iSUB);
            end
        end
    end
end
iEL_beam_end = iEL;
iEL_end = iEL;
fprintf(fileID,'\n');

%% BEAMLOAD    ; Element Beam Loads
fprintf(fileID,'*BEAMLOAD    ; Element Beam Loads\n');
fprintf(fileID,'; ELEM_LIST, CMD, TYPE, DIR, bPROJ, [ECCEN], [VALUE], GROUP\n; ELEM_LIST, CMD, TYPE, TYPE, DIR, VX, VY, VZ, bPROJ, [ECCEN], [VALUE], GROUP\n; [VALUE]       : D1, P1, D2, P2, D3, P3, D4, P4\n; [ECCEN]       : bECCEN, ECCDIR, I-END, J-END, bJ-END\n; [ADDITIONAL]  : bADDITIONAL, ADDITIONAL_I-END, ADDITIONAL_J-END, bADDITIONAL_J-END\n');

ELE_CMD = 'BEAM'; ELE_TYPE = 'UNILOAD'; ELE_DIR = 'GZ'; ELE_NO = 'NO'; ELE_aDir = 'aDir[1]';
ELE_P = '-3e-003';

iEL = iEL_beam_0;
while iEL < iEL_end
    iEL = iEL+1;
    fprintf(fileID,'   %d, %s, %s, %s, %s, %s, %s, , , , 0, %s, 1, %s, 0, 0, 0, 0, , %s, 0, 0, %s, \n',...
        iEL, ELE_CMD, ELE_TYPE, ELE_DIR, ELE_NO, ELE_NO, ELE_aDir,...
        ELE_P, ELE_P, ELE_NO, ELE_NO);
end

%% CONSTRAINT
fprintf(fileID,'*CONSTRAINT    ; Supports\n');
fprintf(fileID,'; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP\n');

iNO = iNO_init; % 初始化iNO
NODE_LIST = sprintf('%dto%d', iNO+1, iNO+XYcoor_num);
CONSTRAINT = 111111; % 6个自由度全约束
fprintf(fileID,'   %s, %d, \n',...
    NODE_LIST, CONSTRAINT);
fprintf(fileID,'\n');

end