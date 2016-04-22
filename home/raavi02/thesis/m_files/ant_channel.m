% ##############################################################################
% ##  ant_channel.m : Simulation von Mobilfunkkanaelen                        ##
% ##############################################################################
%
% Aufruf:       [signal_out,channel] = ant_channel(signal_in,channel);
%
% Version 1.13
%
% Eingabe Parameter signal_in : Einganssignal
%
%   channel             :  structure with all simulation parameter
%
%   channel.name        : Channel Profile Name
%              'h2'        HIPERLAN/2 indoor channels
%              'exp'      Exponential power delay profile/ Classical Doppler 
%              'cost207'  COST207 (GSM)-Profile
%              'lpath'    L-Path Rayleigh Channel 
%                         (all taps are equal distributed)
%     
%
%   channel.number      : h2: for HIPERLAN/2 channel model: 1..7
%                         cost207: and for COST207 power profile: 
%                                    'TU','HT','BU'
%                         lpath:   number of taps (>=1)
%   channel.profile     : 'p' use power profile (default)
%                       : 'd' use delay profile (only 'exp','cost207','lpath')
%   channel.n_echo      : number of echos (default: n_echo = 100)
%   channel.delta_tau   : [s]   power delay spread (only 'exp')
%   channel.tau_max     : [s]   maximum echo delay (only 'exp')
%   channel.f_0         : [Hz]  carrier frequency 
%   channel.v_0         : [m/s] maximum object speed (for doppler spreading)
%                             required for COST207 
%   channel.f_a         : [Hz]  sample rate
%   channel.pause       : [s] pause time, before time variant
%                             convolution (default 0)
%
%  optional:
%   channel.const_taps  : number of non time variant samples (speed up)
%   channel.flag_init   : set if new initialisation is needed
%   channel.seed        : set SEED for rand - function 
% 
%
% output   signal_out
%   channel.echo_state       : current state of each path (vector)
%   channel.echo_doppler     : complex doppler factor of each echo
%   channel.tap_delay        : delay of each tap
%    channel.n_echo          : used number of echos
%   channel.impulse_response : current channel impulse response
%   channel.tau              : tau vector
%   channel.time             : time vector
%   channel.ident            : identification string 
%

% used mex-files: tvc.mex 
%
% ------------------------------------------------------------------------------
%
%  Author:  Heiko Schmidt (University of Bremen)
%  Project:  ---
%  Date:  13-Aug-1999
%  Last modified:   25-Apr-2001 (1.13)
%  Matlab Version: 5.1
%
%  Copyright (c) Department of Telecommunications / University of Bremen 1999
%  This software is property of University of Bremen. Unauthorized 
%  duplication and disclosure to third parties is forbidden. 
% ------------------------------------------------------------------------------
%

function [signal,channel] = ant_channel(signal,channel);
%
signal  = signal(:);
len_sig = length(signal);
%
if isfield(channel,'pause')~=1
  channel.pause = 0;
end;
if (isfield(channel,'flag_init')~=1)
  channel.flag_init = 1;
end;
%
%
% -----------------------------     ----------------------------
%     CHANNEL INITIALIZATION
% -----------------------------     ----------------------------
if channel.flag_init

  % -------------------- SYNTAX CHECK -------------------------
  if isfield(channel,'name')~=1
    channel.name = 'awgn';
  end;
  if isfield(channel,'seed')~=1
    channel.seed = 0;
  end;
  if isfield(channel,'profile')~=1
    channel.profile = 'p';
  end;
  if isfield(channel,'number')~=1
    channel.number = 1;
  end;
  if isfield(channel,'f_0')~=1
    channel.f_0   = 5.2e9;
  end;
  if isfield(channel,'v_0')~=1
    channel.v_0   = 0;
  end;
  if isfield(channel,'f_a')~=1
    channel.f_a   = 1;
  end;
  if isfield(channel,'n_echo')~=1
    channel.n_echo = 300;
  end;
  if isfield(channel,'const_taps')~=1
    channel.const_taps = 1;
  end;
  if isfield(channel,'flag_show')~=1
    channel.flag_show = 0;
  end;

  if isfield(channel,'tau_max')~=1
    channel.tau_max = 1;
  end;

  if isfield(channel,'delta_tau')~=1
    channel.delta_tau = 1;
  end;

  if isfield(channel,'pause')~=1
    channel.pause = 0;
  end;


  % ----------------- END OF SYNTAX CHECK ---------------------
  c_0 = 3e8;
  if channel.seed
    rand('seed',channel.seed);
  end;
  channel.memory = 0;
  %
  %
  % -----------------------------------------------------------
  switch upper(channel.name)
  case 'H2'
    % --------------- ETSI EP BRAN - HIPERLAN/2 INDOOR CHANNEL MODELS ----------
    switch channel.number 
    case 1 % channel modell A (50 ns)
      v_arp    = 10.^([0.0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 ...
                       -4.7 -7.3 -9.9 -12.5 -13.7 -18.0 -22.4 -26.7]/20);
      v_delay  = [0 10 20 30 40 50 60 70 80 90 110 140 170 200 240 290 ...
                  340 390] * 1e-9;
      rice=0;
    case 2 % channel modell B (100 ns)
      v_arp    = 10.^([-2.6 -3.0 -3.5 -3.9 -0.0 -1.3 -2.6 -3.9 -3.4 -5.6 ...
                       -7.7 -9.9 -12.1 -14.3 -15.4 -18.4 -20.7 -24.6]/20);
      v_delay  = [0 10 20 30 50 80 110 140 180 230 280 330 380 430 490 ...
                  560 640 730] * 1e-9;
      rice=0;

    case 3  % channel modell C (150 ns)
      v_arp    = 10.^([-3.3 -3.6 -3.9 -4.2 -0.0 -0.9 -1.7 -2.6 -1.5 -3.0 ...
                       -4.4 -5.9 -5.3 -7.9 -9.4 -13.2 -16.3 -21.2]/20);
      v_delay  = [0 10 20 30 50 80 110 140 180 230 280 330 400 490 600 ...
                  730 880 1050] * 1e-9;
      rice=0;

    case 4  % channel modell D (150 ns)
      v_arp    = 10.^([-10.0 -10.3 -10.6 -6.4 -7.2 -8.1 -9.0 -7.9 -9.4 ...
                       -10.8 -12.3 -11.7 -14.3 -15.8 -19.6 -22.7 -27.6]/20);
      v_delay  = [10 20 30 50 80 110 140 180 230 280 330 400 490 600 ...
                  730 880 1050] * 1e-9;
      rice     = 0.75^2;

    case 5  % channel modell E (250 ns)
      v_arp    = 10.^([-4.9 -5.1 -5.2 -0.8 -1.3 -1.9 -0.3 -1.2 -2.1 0.0 ...
                       -1.9 -2.8 -5.4 -7.3 -10.6 -13.4 -17.4 -20.9]/20);
      v_delay  = [0 10 20 40 70 100 140 190 240 320 430 560 710 880 ...
                  1070 1280 1510 1760] * 1e-9;
      rice=0 ; % Linear !!!;

    case 6  % channel modell F (10 ns)
      v_arp    = 10.^([0.0  -4.3  -8.7 -13.0 -17.4 -21.7 -26.1 -30.4 -34.7 ...
                       -39.1 -43.4 -47.8 -52.1 -56.5 -60.8 -65.1 -69.5 ...
                       -73.8 ]/20); 
      v_delay = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 ...
                 170 ] *1e-9; 
      rice=0;

    case 7  % channel modell G (30 ns)
      v_arp   = 10.^([ 0.0  -1.4  -2.9  -4.3  -5.8  -7.2  -8.7 -10.1 ...
                      -11.6 -13.0 -14.5 -13.6 -16.5 -19.4 -22.3 -25.2 ...
                      -28.1 -29.8 ]/20); 
      v_delay = [0 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200 220 ...
                 250 ] *1e-9; 
      rice=0;
    otherwise 
      v_arp  = 1;
      v_delay  = 0;
      rice     = 0;
    end; 

    tap_delay_tmp = round(v_delay * channel.f_a);
    n_tap      = length(v_delay);
    n_echo_t      = ceil(channel.n_echo/n_tap);
    channel.n_echo= n_echo_t * n_tap;
    total_power  = sum(abs(v_arp).^2);
    norm_power  = 1 / sqrt(n_echo_t * total_power);
    w_doppler  = 2 * pi * channel.v_0 / c_0 * channel.f_0 / channel.f_a;
    echo_state_tmp        = exp(i * 2 * pi * rand(n_echo_t,n_tap)) ...
      * diag(v_arp) * norm_power;
    echo_doppler_tmp      = exp(i * w_doppler* ...
                                cos(2 * pi * rand(n_echo_t,n_tap) ));
    tap_delay_tmp  = repmat(tap_delay_tmp,n_echo_t,1);
    channel.tap_delay  = tap_delay_tmp(:);
    channel.echo_state  = echo_state_tmp(:);
    channel.echo_doppler  = echo_doppler_tmp(:);
    if rice > 0
      %rice_lin  = 10^(rice/20);
      channel.echo_state   = [sqrt(rice);channel.echo_state];
      channel.echo_doppler = [1;channel.echo_doppler];
      channel.state(1)     = exp(i*2*pi*rand(1,1))*sqrt(rice);
      norm          = 1/sqrt(1+rice);
      channel.echo_state   = channel.echo_state * norm;
      channel.tap_delay    = [0;channel.tap_delay];
    end; % END RICE
    channel.ident = sprintf('H2n%dne%dv%dfa%d',...
                            channel.number,round(channel.n_echo/10),...
                            round(channel.v_0*10),round(channel.f_a/1000000));
    % ------------------------------ END OF HIPERLAN/2 -------------------------
    %
    %
    % --------------------------------------------------------------------------
    %
    % --------------------------------- COST 207 PROFILES ----------------------
    %
    % --------------------------------------------------------------------------
  case 'COST207'
    %
    switch upper(channel.profile)
      %
      % ---------- COST 207 PROFILES with power profiles -----------------------
      %
    case 'P'
      %
      %
      switch upper(channel.number)
      case 'TU'
        channel.tau_max   = 7e-6;
        channel.delta_tau = 1e-6;
        n_tap            = ceil(channel.tau_max * channel.f_a);
        n_echo_t        = ceil(channel.n_echo/n_tap);
        channel.n_echo        = n_echo_t * n_tap;
        tau            = [0:n_tap-1] / channel.f_a;
        power_profile         = sqrt(exp(-tau/channel.delta_tau));
        channel.ident = sprintf('TU%dpp',round(3.6*channel.v_0));
      case 'BU' 
        channel.tau_max   = 10e-6;
        channel.tau_max1  = 5e-6;
        channel.delta_tau = 1e-6;
        n_tap    = ceil(channel.tau_max * channel.f_a);
        n_tap1            = ceil(channel.tau_max1 * channel.f_a);
        n_tap2    = n_tap - n_tap1;   
        n_echo_t          = ceil(channel.n_echo/n_tap);
        channel.n_echo    = n_echo_t * n_tap;
        tau1           = [0:(n_tap1-1)] / channel.f_a;
        tau2           = [n_tap1:(n_tap-1)] / channel.f_a;
        power_profile1    = sqrt(exp(-tau1/channel.delta_tau));
        power_profile2    = sqrt(0.5 * exp(-(tau2-channel.tau_max1)...
                                           /channel.delta_tau));
        tau            = [tau1 tau2];
        power_profile     = [power_profile1 power_profile2];
        channel.ident = sprintf('BU%dpp',round(3.6*channel.v_0));
      case 'HT'
        channel.tau_max = 20e-6;
        channel.tau_max1= 2e-6; 
        channel.tau_max2= 15e-6;
        channel.delta_tau = 1e-6;
        n_tap1          = ceil(channel.tau_max1 * channel.f_a);
        n_tap2  = ceil((channel.tau_max-channel.tau_max2) * channel.f_a);
        n_tap  = n_tap1 + n_tap2;
        n_echo_t      = ceil(channel.n_echo/n_tap);
        channel.n_echo  = n_echo_t * n_tap;
        tau1  = [0:(n_tap1-1)] / channel.f_a;
        tau2  = channel.tau_max2 + [0:(n_tap2-1)] / channel.f_a;
        power_profile1  = sqrt(exp(-tau1/channel.delta_tau));
        power_profile2  = sqrt(0.04 * exp(-(tau2-channel.tau_max2) ...
                                          /channel.delta_tau));
        tau   = [tau1 tau2];
        power_profile   = [power_profile1 power_profile2];
        channel.ident   = sprintf('HT%dpp',round(3.6*channel.v_0));
      end;
      total_power       = sum(abs(power_profile).^2);
      norm_power        = 1 / sqrt(n_echo_t * total_power);
      echo_state_tmp    = exp(i * 2 * pi * rand(n_echo_t,n_tap)) ...
        * diag(power_profile) * norm_power;
      tap_delay_tmp     = repmat([round(tau*channel.f_a)],n_echo_t,1);
      channel.tap_delay = tap_delay_tmp(:);
      channel.echo_state= echo_state_tmp(:);
      %
      %
      % ------------------------ END OF POWER PROFILE -------------------------
      %
      % --------------------------- POWER DELAY PROFILE ------------------------
      %
    case 'D'
      %
      %
      switch upper(channel.number)
      case 'TU'
        channel.tau_max    = 7e-6;
        channel.delta_tau  = 1e-6;
        max_limit    = exp(-channel.tau_max/channel.delta_tau);
        n_tap      = channel.n_echo;
        equal_dist       = (1-max_limit) * rand(n_tap,1) + max_limit;
        delays   = round(-log(equal_dist) * (channel.f_a*channel.delta_tau));
        channel.ident     = sprintf('TU%ddp',round(3.6*channel.v_0));
      case 'BU'
        channel.tau_max   = 10e-6;
        channel.tau_max1  = 5e-6;
        channel.delta_tau = 1e-6;
        min_limit2 = 0.5;
        max_limit1 = exp(-channel.tau_max1/channel.delta_tau);
        max_limit2 = exp(-(channel.tau_max-channel.tau_max1)...
                         /channel.delta_tau) * min_limit2;
        n_tap      = channel.n_echo;
        dist_factor= (1-max_limit1) + (min_limit2*(1-max_limit2));
        equal_dist = dist_factor * rand(n_tap,1) -(min_limit2*(1-max_limit2));
        index1     = find(equal_dist > 0);
        index2     = find(equal_dist < 0);
        delays1    = round(-log(equal_dist(index1)+max_limit1) ...
                           * (channel.f_a*channel.delta_tau));
        delays2    = round(-log((max_limit2 ...
                                 - equal_dist(index2))/min_limit2) ...
                           * (channel.f_a*channel.delta_tau) ...
                           + (channel.f_a*channel.tau_max1));
        delays     = [delays1;delays2];                
        channel.ident     = sprintf('BU%ddp',round(3.6*channel.v_0));
      case 'HT'
        channel.tau_max   = 20e-6;
        channel.tau_max1  = 2e-6;
        channel.tau_max2  = 15e-6;

        channel.delta_tau = 1e-6;
        min_limit2    = 0.04;
        max_limit1    = exp(-channel.tau_max1/channel.delta_tau);
        max_limit2    = exp(-(channel.tau_max-channel.tau_max2) ...
                            /channel.delta_tau) * min_limit2;
        n_tap      = channel.n_echo;
        dist_factor= (1-max_limit1) + (min_limit2*(1-max_limit2));
        equal_dist = dist_factor * rand(n_tap,1) -(min_limit2*(1-max_limit2));
        index1     = find(equal_dist > 0);
        index2     = find(equal_dist < 0);
        delays1       = round(-log(equal_dist(index1)+max_limit1) ...
                              * (channel.f_a*channel.delta_tau));
        delays2       = round(-log((max_limit2 ...
                                    - equal_dist(index2))/min_limit2) ...
                              * (channel.f_a*channel.delta_tau) ...
                              + (channel.f_a*channel.tau_max2));
        delays     = [delays1;delays2];        
        channel.ident   = sprintf('HT%ddp',round(3.6*channel.v_0));
      end;
      channel.tap_delay  = delays;
      channel.echo_state = exp(i * 2 * pi * rand(n_tap,1)) /sqrt(n_tap); 
      %
      % ------------------------ END OF POWER DELAY PROFILE --------------------
    end;
    %
    % ----------------- compute doppler frequencies --------------------
    %
    n_doppler  = length(channel.tap_delay);

    limit1       = 500e-9 * channel.f_a ;
    limit2          = 2e-6   * channel.f_a;
    f_doppler = channel.v_0 / c_0 * channel.f_0;
    f111 = -0.8 * f_doppler; % tau < 2 us
    f112 = 0.05 * f_doppler; % tau < 2 us
    f121 = 0.4  * f_doppler; % tau < 2 us
    f122 = 0.1  * f_doppler; % tau < 2 us

    f211 = 0.7  * f_doppler; % tau > 2 us
    f212 = 0.1  * f_doppler; % tau > 2 us
    f221 = -0.4 * f_doppler; % tau > 2 us
    f222 = 0.15 * f_doppler; % tau > 2 us

    index_class = find(channel.tap_delay <= limit1);
    index_gaus1 = find(channel.tap_delay > limit1 ...
                       & channel.tap_delay <= limit2);
    index_gaus2 = find(channel.tap_delay > limit2);

    n_class = length(index_class);
    n_gaus1 = length(index_gaus1);
    n_gaus2 = length(index_gaus2);

    % ---- classical doppler spectrum 
    class_doppler = exp((2 * pi * i * f_doppler/channel.f_a) ...
                        * cos(2 * pi * rand(n_class,1)));
    %
    %
    % ---- Gauss doppler spectrum tau < 2 us
    %
    f_d_tmp =  [(f112)*randn(n_gaus1,1)+(f111) ...
                (f122)*randn(n_gaus1,1)+(f121)];    
    %
    if f_doppler ~=0
      v_dec   =  ([1:n_gaus1]')  + (n_gaus1 ...
      * sign(round(1.055 * sqrt(f122/f112) * rand(n_gaus1,1))));
      f_d_tmp = f_d_tmp(v_dec);
    else 
      f_d_tmp = f_d_tmp(1:n_gaus1);    
    end;
    %
    gaus1_doppler = exp((2 * pi * i/channel.f_a) * f_d_tmp);
    %
    %
    % ---- Gauss doppler spectrum tau > 2 us
    %
    f_d_tmp =  [(f212)*randn(n_gaus2,1)+(f211) ...
                (f222)*randn(n_gaus2,1)+(f221)];    
    %
    if (f_doppler ~=0)
      v_dec   =  ([1:n_gaus2]') + (n_gaus2 ...
      * sign(round(1.03 * sqrt(f222/f212)...
                   * rand(n_gaus2,1))));
      f_d_tmp = f_d_tmp(v_dec);
    else 
      f_d_tmp = f_d_tmp(1:n_gaus2);    
    end;
    %        
    gaus2_doppler = exp((2 * pi * i/channel.f_a) * f_d_tmp);
    %
    % ------------------------------------------------
    %
    channel.echo_doppler = zeros(n_doppler,1);
    channel.echo_doppler(index_class) = class_doppler;
    channel.echo_doppler(index_gaus1) = gaus1_doppler;
    channel.echo_doppler(index_gaus2) = gaus2_doppler;
    %
    %
    %
    % ------------------------ END od COST 207 -----------------------
    %
    %
    % -------------------------- L-Path Channel Model --------------------------
    %
  case 'LPATH'
    switch upper(channel.profile)
    case 'P' 
      w_doppler  = 2 * pi * channel.v_0 / c_0 * channel.f_0 / channel.f_a;
      n_tap   = channel.number;
      n_echo_t      = ceil(channel.n_echo/n_tap);
      channel.n_echo  = n_echo_t * n_tap;
      tau  = [0:n_tap-1] / channel.f_a;
      norm_power  = 1 / sqrt(n_echo_t * n_tap);
      echo_state_tmp     = exp(i * 2 * pi * rand(n_echo_t,n_tap))  * norm_power;
      echo_doppler_tmp   = exp(i * w_doppler* ...
                               cos(2 * pi * rand(n_echo_t,n_tap) ));
      tap_delay_tmp        = repmat([0:n_tap-1],n_echo_t,1);

      channel.tap_delay  = tap_delay_tmp(:);
      channel.echo_state  = echo_state_tmp(:);
      channel.echo_doppler   = echo_doppler_tmp(:);

      channel.ident = sprintf('LP%d',channel.number);
      %
      % --------- Equal DISTRIBUTED POWER DELAY PROFILE ------------------------
      %
    case 'D'
      max_limit    = channel.number;
      n_echo      = channel.n_echo;
      w_doppler     = 2 * pi * channel.v_0 / c_0 * channel.f_0 / channel.f_a;
      equal_dist       = floor(max_limit * rand(n_echo,1));
      channel.tap_delay  = equal_dist;
      channel.echo_state = exp(i * 2 * pi * rand(n_echo,1)) /sqrt(n_echo);
      channel.echo_doppler = exp(i * w_doppler* cos(2 * pi * rand(n_echo,1) ));

      channel.ident = sprintf('LD%d',channel.number);
    end;

    %
    % -------------------------- End of L-Path --------------------------
    %
    % ---------- EXPONENTIALLY DISTRIBUTED POWER PROFILE -----------------------
    %
  case 'EXP'
    switch upper(channel.profile)
    case 'P' 
      w_doppler  = 2 * pi * channel.v_0 / c_0 * channel.f_0 / channel.f_a;
      n_tap   = ceil(channel.tau_max * channel.f_a);
      n_echo_t      = ceil(channel.n_echo/n_tap);
      channel.n_echo  = n_echo_t * n_tap;
      tau  = [0:n_tap-1] / channel.f_a;
      power_profile   = sqrt(exp(-tau/channel.delta_tau));
      total_power  = sum(abs(power_profile).^2);
      norm_power  = 1 / sqrt(n_echo_t * total_power);
      echo_state_tmp        = exp(i * 2 * pi * rand(n_echo_t,n_tap)) ...
        * diag(power_profile) * norm_power;
      echo_doppler_tmp      = exp(i * w_doppler ...
                                  * cos(2 * pi * rand(n_echo_t,n_tap) ));
      tap_delay_tmp  = repmat([0:n_tap-1],n_echo_t,1);

      channel.tap_delay  = tap_delay_tmp(:);
      channel.echo_state  = echo_state_tmp(:);
      channel.echo_doppler   = echo_doppler_tmp(:);

      channel.ident = sprintf('EPdt%d',round(1e9*channel.delta_tau));
      %
      % ----- EXPONENTIALLY DISTRIBUTED POWER DELAY PROFILE --------------------
      %
    case 'D'
      max_limit    = exp(-channel.tau_max/channel.delta_tau);
      n_tap      = channel.n_echo;
      w_doppler     = 2 * pi * channel.v_0 / c_0 * channel.f_0 / channel.f_a;
      equal_dist       = (1-max_limit) * rand(n_tap,1) + max_limit;
      delays       = round(-log(equal_dist) * (channel.f_a*channel.delta_tau));
      channel.tap_delay  = delays;
      channel.echo_state = exp(i * 2 * pi * rand(n_tap,1)) /sqrt(n_tap);
      channel.echo_doppler = exp(i * w_doppler* cos(2 * pi * rand(n_tap,1) ));

      channel.ident = sprintf('EDdt%d',round(1e9*channel.delta_tau));
    end;
    %
    % --------------------------------------------------------------------------
    %
  case 'TEST'

    channel.ident = 'TesT';
    channel.tap_delay = [0 0 1 1 2]';

    channel.echo_doppler   = exp(2*pi*j*[0 0.0015 0 -0.0015 0]');
    channel.echo_state     = [1 0.9 0.6 j*0.5 j*0.5]';
    channel.n_echo         = 5;
    channel.const_taps     = 1;
    channel.v_0 = 100;

    % --------------------------------------------------------------------------
    %

  otherwise
    channel.ident = '';
    channel.echo_doppler = 1;
    channel.echo_state   = 1;
    channel.v_0          = 0;
    channel.tap_delay    = 0;
  end;
  if channel.seed > 0 
    channel.ident = sprintf('%sSEED%d',channel.ident,channel.seed);
  end;

  % ---------------- IMPORTANT for all channel types ! ----------------------
  if channel.const_taps > 1
    channel.echo_doppler = channel.echo_doppler .^ channel.const_taps;
  elseif (channel.const_taps == 0) | (channel.v_0 == 0);
    channel.echo_doppler = ones(size(channel.echo_doppler));
    channel.const_taps = 0;
    channel.v_0 = 0;
  end;
  if (channel.v_0 == 0)
    channel.const_taps = 0;
  end;

  channel.flag_init = 0;
end;
% -----------------------------     ----------------------------
%     END OF CHANNEL INITIALIZATION
% -----------------------------     ----------------------------
%
% --------------------------------------------------------------
% ------------------------ PAUSE SIMULATION --------------------
% --------------------------------------------------------------
%
if (channel.pause > 0)
  n_pause = round(channel.pause * channel.f_a);
  channel.echo_state = channel.echo_state .* (channel.echo_doppler ...
                                              .^ (n_pause/channel.const_taps));
  mem_diff = length(channel.memory) - n_pause;
  if (mem_diff <= 0)
    channel.memory = 0;
  else
    channel.memory = [zeros(mem_diff,1) ;...
                      channel.memory(n_pause+1:n_pause+mem_diff)];
  end;
  channel.pause = 0;
end;
%
% ------------------- END OF PAUSE ---------------
%
echo_state_show = channel.echo_state;
% 
% ------------------------------------------------------------------------------
% -------------------------    time variant convolution   ----------------------
% ------------------------------------------------------------------------------

tap_delay  = channel.tap_delay;
offset     = max(channel.tap_delay);
sig_len    = length(signal);
%
if sig_len < (offset+1)
  signal  = [signal(:);zeros(offset-sig_len+1,1)];
  sig_len = length(signal);
end;
%
if (channel.memory(1) ~= 0)
  mem_len    = length(channel.memory);
  signal_in  = [[zeros(max(offset-mem_len,0),1);...
                 channel.memory(1:min(offset,mem_len))];signal];
else
  signal_in  = [zeros(offset,1);signal];
end;
channel.memory  = signal(sig_len-offset+1:sig_len);  % save input signal
%                                                      for next convolution 
v_state     = channel.echo_state(:);
v_doppler   = channel.echo_doppler(:);
[n_tap]     = length(channel.echo_state);

%
%
%
[signal_out,v_state_back,v_ipr] = tvc(signal_in,offset,v_state,...
                                      v_doppler,channel.const_taps,...
                                      n_tap,channel.tap_delay);
%
%
%
l_ipr    = offset + 1;
n_ipr    = floor(length(v_ipr) / l_ipr);
m_ipr    = zeros(l_ipr,n_ipr);
m_ipr(:) = v_ipr(1:(l_ipr*n_ipr));

% ---------------- COMPUTE TIME VARIANT IMPULSE RESPONSE -----------------------
%
% ------------------------------------------------------------------------------
channel.impulse_response = m_ipr;
channel.tau  = ([0:l_ipr-1].') * (1/channel.f_a);
channel.time = ([0:n_ipr-1]) * (1/channel.f_a * channel.const_taps);

signal = signal_out;
channel.echo_state = v_state_back;

% ### EOF ######################################################################



