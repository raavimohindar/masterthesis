% Function file: Symbol by Symbol MAX-LOG-MAP decoding of Walsh-Chips
%% Copyright 1998 Armin Dekorsy
% for academic use only
%
%function [L_out]=max_walsh(L_xy,L_in,st_walsh,m_ary,n_sym)
%
% -------------------------------------------------------------------------------------------------
% PROJECT:   Channel-Coding-Algorithmn
% Submitted by: University of Bremen
%  
% -------------------------------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------------------------------
% 1. Symbol by Symbol Max-Log-MAP-Decoding of Walsh-Chips (Walsh-Symbols of Walsh-Codewords)
% 2. Direct Implementation of MAP-Criteria by using DHT
% 3. Delivers Log-Likelihood-Ratios (LLR) of systematic symbols (infobits) due to MAP-criteria
% 4. Input-Signals are LLR`s 
% -------------------------------------------------------------------------------------------------
% INPUT:
% -------------------------------------------------------------------------------------------------
% L_xy:     L_xy=LLR(y/x)=4*snr*y=L_c*y with L_c equals Channel-State-Information
%           m_ary*#of transmitted symbols (codewords)
% L_in:     A-priori LLR of systematic symbols (infobits)
%           ldM * # of transmitted symbols (codewords)
% st_walsh: structure of walsh-informations: see [st_walsh]=walsh_par(m_ary)
% m_ary:    M-ary orthogonal Modulation
% n_sym:    # of transmitted symbols (codewords)
% -------------------------------------------------------------------------------------------------
% OUTPUT:
% -------------------------------------------------------------------------------------------------
% - L_out: LLR of systematic symbols (infobits), ldM*# of transmitted (symbols) codewords  
% -------------------------------------------------------------------------------------------------
% FUNCTIONS:
% -------------------------------------------------------------------------------------------------
% - to obtain st_walsh use walsh_par.m
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
% Author: Armin Dekorsy  17.12.98 
%% -------------------------------------------------------------------------------------------------

function [L_out]=max_walsh(L_xy,L_in,st_walsh,m_ary,n_sym)

ldM=log2(m_ary);
L_a=zeros(m_ary,n_sym);
L_out=zeros(ldM,n_sym);

L_a(st_walsh.sys_k,:)=L_in;

LLR=L_xy+L_a;

DHT_llr = st_walsh.W*LLR;             %Berechnung der A-posteriori-Symbolw'keiten
                                         %durch LLR(Y/X) ausgedrueckt 
W_part=st_walsh.W(st_walsh.sys_k,:);

for run0=1:ldM
  l=find(W_part(run0,:)==1);
  Max1=max(DHT_llr(l,:));
  l=find(W_part(run0,:)==-1);
  Max_1=max(DHT_llr(l,:));
  L_out(run0,:)=0.5*(Max1-Max_1);       %inneren Infobit (systematischen Bits): 
                                         %LLR=Lc*y+L(u_k)+L_ext
                                         %ldM*n_sym-matrix
end; 

 
                                         
