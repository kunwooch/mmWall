clear;clc;
close all
addpath './functions/'
%% import Huygen's pattern
load("./VNA_mmWall.mat")

freq                            = 24.3;
freq_idx                        = find(freq_range==freq);
S11                             = S11(freq_idx, :, :);
S21                             = S21(freq_idx, :, :);

%% Redefine board for testing 
num_boards = 76;
boards = ones(1, num_boards);

plot_huygens(squeeze(S11(1, :, :)), squeeze(S21(1, :, :)), num_UE, num_UM, freq);

%% Define parameters
[param, eval] = define_params(freq);
disp(strcat("Gain limit: ", num2str(param.gain_dB)))

param.half                      = 0;
param.mode                      = 'transmission';     % transmission or reflection or bidirection or combine
param.dual                      = 0;                  % transmission or reflection or bidirection or combine
param.num_UE                    = num_UE;
param.num_UM                    = num_UM; 
param.boards                    = boards; 
param.S11                       = S11;
param.S21                       = S21;

%% Change in steering angle
Rx_angle1                   = [-70:35:70];
Tx_angle1                   = 0;

for angle_idx = 1:length(Rx_angle1)  
     %% Single-beam
    fprintf("Tx Angle: %d, Rx Angle: %d \n", Tx_angle1, Rx_angle1(angle_idx))

    [C_gt, ue, um, phi_nm, ~]          = optimize_mmWall(Tx_angle1, Rx_angle1(angle_idx), param, eval);

    [C_tra, err]                       = volt2coef(ue, um, S11, S21, num_UE, num_UM, boards, 'transmission');
    [C_ref, err]                       = volt2coef(ue, um, S11, S21, num_UE, num_UM, boards, 'reflection');

    C_gt                                = repmat(C_gt, [param.n, 1]);
    C_tra                              = repmat(C_tra, [param.n, 1]);
    C_ref                              = repmat(C_ref, [param.n, 1]);

    G_gt                                = abs(sum(sum(C_gt.*exp(1i*phi_nm))))^2;
    if strcmp(param.mode, 'transmission')
        G_opt                           = abs(sum(sum(C_tra.*exp(1i*phi_nm))))^2; 
        C_gt1                           = C_gt;
        C_gt2                           = zeros(size(C_gt));
    else
        G_opt                           = abs(sum(sum(C_ref.*exp(1i*phi_nm))))^2;
        C_gt1                           = zeros(size(C_gt));
        C_gt2                           = C_gt;
    end

    fprintf("Perfect gain: %d, mmWall gain (opt): %d\n", pow2db(G_gt)/2, pow2db(G_opt)/2)
    [af_gt, af_opt] = array_factor(C_gt1, C_gt2, C_tra, C_ref, param, 0);

    graph_beamsteering(Tx_angle1, Rx_angle1(angle_idx), C_tra, C_ref, C_gt1, C_gt2, af_opt, af_gt, param);
end
