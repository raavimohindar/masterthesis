% Function file: Calculation of SS-MAP-Decoding parameters of Walsh-Codes
%% Copyright 1998 Armin Dekorsy
% for academic use only
%
% function [st_walsh]=walsh_par(m_ary)
%
% -------------------------------------------------------------------------------------------------
% PROJECT:   Channel-Coding-Algorithmn
% Submitted by: University of Bremen
%  
% -------------------------------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------------------------------
% 1. Calculates necessary parameters for SS-MAP-Decoding of Walsh-Codes
% 2. Used in function ss_map_walsh.m 
% -------------------------------------------------------------------------------------------------
% INPUT:
% -------------------------------------------------------------------------------------------------
% m_ary:    M-ary orthogonal Modulation
% -------------------------------------------------------------------------------------------------
% OUTPUT:
% -------------------------------------------------------------------------------------------------
% - st_walsh: Structur-Array with following elements:
%             W:        m_ary Hadamard-Matrix
%             W_part1:  part of Hadamard-Matrix where systematic symbols ==1
%                       ldM*m_ary-matrix
%             W_part_1: part of Hadamard-Matrix where systematic symbols ==-1
%                       ldM*m_ary-matrix
%             sys_k: Positions of systemtic symbols in codeword  
% -------------------------------------------------------------------------------------------------
% FUNCTIONS:
% -------------------------------------------------------------------------------------------------
% - used in ss_map_walsh.m
% -------------------------------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------------------------------
% - see Iterative Decoding and Despreading improves CDMA-Systems using
%   M-ary Orthogonal Modulation and FEC
%   R. Herzog, A. Schmidbauer, J. Hagenauer
%   ICC-97, pp.909-913
%   Equation (2)
%      
% -------------------------------------------------------------------------------------------------
% IMPORTANT PARAMETERS:
% -------------------------------------------------------------------------------------------------
% - none
%% -------------------------------------------------------------------------------------------------
% Author: Armin Dekorsy  04.12.98 
%% -------------------------------------------------------------------------------------------------


function [st_walsh]=walsh_par(m_ary)

ldM=log2(m_ary);
st_walsh.W=hadamard(m_ary);               %Hadamard-Matrix
k=ldM:-1:1;                               %Addition der zu einem systematischen Bit 
st_walsh.sys_k=(m_ary./(2.^k))+1;         %zugehoerigen A-posteriori-Symbolw`keiten
W_part=st_walsh.W(st_walsh.sys_k,:);      %ldM*m_ary-Matrix with inner systematic bit
st_walsh.W_part1=(W_part==1);             %ldM*m_ary Matrix with ones where inner systematic bit=1,
                                          %otherwise zero
st_walsh.W_part_1=(W_part==-1);           %ldM*m_ary Matrix with one  where inner systematic bit=-1, 
                                          %otherwise zero
