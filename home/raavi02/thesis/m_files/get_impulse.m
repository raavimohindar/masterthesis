% ##############################################################################
% ##  get_impulse : Erzeugung von beliebig abgetasteten Impulsantworten       ##
% ##############################################################################
%
% Erzeugung von beliebig abgetasteten Impulsantworten von Cosinus-Roll-
% off-, Wurzel-Cosinus-Roll-off-, Gauss- und GMSK-(Gauss*Rect)-Filtern
% (Zeitbereichs-Entwuerfe).
%
% function [gn, tn, imp_name] = ...
%             get_impulse (imp_type, imp_par, M_over, LT[, a[, norm_f]])
% ----------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------
% imp_type/  String, aus dem sich die gewuenschte Impulsform und die
% imp_par :  Bedeutung von imp_par ergibt:
%            -----------------------------------------------------------
%            | imp_type |          Impulstyp                 | imp_par |
%            | enthaelt |                                    | bedeutet|
%            -----------------------------------------------------------
%            | 'SQRT'   | Wurzel-Cosinus-Roll-off-Impuls     | rolloff |
%            | 'WUR'    |           -- " --                  | rolloff |
%            | 'COS'    | Cosinus-Roll-off-Impuls            | rolloff |
%            | 'GMSK'   | Gauss-Impuls gefaltet mit rect(t/T)| f3dB_T  |
%            | 'RECT'   |           -- " -- (!)              | f3dB_T  |
%            | 'GAU'    | Gauss-Impuls                       | f3dB_T  |
%            -----------------------------------------------------------
%            wobei fuer - den Roll-off Faktor gilt  : 0 <= rolloff <= 1,
%                       - die norm. 3dB-Bandbr. gilt: f3dB_T:=f3dB*T>=0.
% M_over  :  Anzahl der Abtastwerte pro Symbolperiode T.
% LT      :  Gewuenschte einseitige(!) Impulslaenge in Symbolperioden T.
% [a]     :  Empirischer Faktor zur Skalierung der Zeitachse des Impulses
%            mit der Wirkung, dass die Nyquistfreq. des Filters ein wenig
%            verkleinert wird. Nach der anschliessenden Fensterung soll
%            sie dann wieder den korrekten Wert haben.
%            Defaultwert: a=1.0 => keine Skalierung der Zeitachse.
% [norm_f]:  Optionales flag, das die Normierung des outputs gn angibt:
%      ~=1:  Normiere im Zeitbereich auf Maximalwert eins (default)
%      ==1:  Norm. im Freq.bereich, so dass max{|Gn|}==1 mit Gn=FFT{gn}.
%            N.B.:  der exakte Wert eins ergibt sich nur fuer LT=inf.
% ----------------------------------------------------------------------
% OUTPUT:
% ----------------------------------------------------------------------
% gn      :  Abtastung des normierten kont. Impulses zu den Zeitpunkten
%            t = (a*n)*Ta = (a*n)*T/M_over mit n = -LT*M_over:LT*M_over.
%            Der Spaltenvektor gn enthaelt die L := 2*LT*M_over+1 so
%            entstandenen Abtastwerte, d.h.
%            - der Hauptkoeffizient liegt auf dem Index LT*M_over+1;
%            - die ersten (letzten) LT*M_over Koeffizienten gehoeren zu
%              negativen (positiven) Zeitpunkten;
%            - Gewinnt man die Impulsantwort im Symboltakt durch
%              gn(1:M_over:2*LT*M_over+1), so hat sie die Laenge 2*LT+1.
% tn      :  Spaltenvektor der Laenge L mit den norm. Abtastzeitpunkten,
%            die zu jedem Wert von gn gehoeren: t = (a*n)*T/M_over =>
%            tn := t/T =a*n/M_over. Mit plot(tn,gn) kann der Impuls als
%            Funktion von t/T dargestellt werden.
% imp_name:  String, der den Impulstyp und den Impulsparameter angibt,
%            z.B. 'SqrtCosRoll (r=0.5)' oder 'GMSK (f3T=%3.1f)'.
% -----------------------------------------------------------------------
% REMARKS:
% -----------------------------------------------------------------------
% - Herleitung Wurzel-Cos-Roll-off:  Wolfgang Behm: "Eine Datenempfaenger
%   struktur mit adaptivem rekursivem Entzerrer und Viterbi-Detektor",
%   Dissertation an der TU Hamburg-Harburg, S. 129ff, Mai 1991.
% - Herleitung Gauss und GMSK:  Kammeyer: "Nachrichtenuebertragung"
% ----------------------------------------------------------------------
% AUTHORS:  Marcus Benthin,     01-jul-92  (Wurzel-Cosinus-Roll-off)
%           Dieter Boss,        14-aug-96  (" und Cosinus-Roll-off)
%           Thorsten Petermann, 12-jul-96  (Gauss und GMSK)
%           Dirk Nikolai,       08-nov-97   ("end" at end of function
%                                            removed for version >= 5.0)
% ----------------------------------------------------------------------

function [gn, tn, imp_name] = get_impulse (imp_type, imp_par,...
                                           M_over, LT, a, norm_f)

% #####  S0:  Impulstyp feststellen
imp_type = imp_type(find(isletter(imp_type)));    % Impulsart
imp_type = upper(imp_type);

if     isempty(findstr(imp_type,'SQRT')) == 0, imp_no = 1;% Erkannte Impulsarten
elseif isempty(findstr(imp_type,'WUR'))  == 0, imp_no = 1;
elseif isempty(findstr(imp_type,'COS'))  == 0, imp_no = 2;
elseif isempty(findstr(imp_type,'GMSK')) == 0, imp_no = 3;
elseif isempty(findstr(imp_type,'RECT')) == 0, imp_no = 3;
elseif isempty(findstr(imp_type,'GAU'))  == 0, imp_no = 4;
else,  error (sprintf('ERROR (get_impulse):  Unknown impulse type "%s"!',...
                      imp_type));
end;

% #####  S1:  Positive norm. Abtastzeitpunkte berechnen
if nargin<5, a = 1.0;  end;               % Default: kein Dehnungsfaktor fuer 
%                                           Zeitachse
if nargin<6, norm_f = 0; end;             % Default: normalize in time domain
%                                           to unit maximum value
n_pos = 1:LT*M_over;                      % Positive diskrete Zeitachse
an_pos = a*n_pos.';                       % Dehnungsfaktor fuer die Zeitachse 
%                                           --> Verlegung des Sym.punktes der
%                                           resultierenden Nyquistflanke
tn_pos = an_pos/M_over;                   % Pos. norm. Abtastzeitpunkte
%                                           t/T = an_pos/M_over.
% tn_pos = \omega_N t/pi fuer t=an_pos T/M_over.


% #####  S2:  Impulskoeff. gn_pos fuer t>0 berechnen
if (imp_no==1) | (imp_no==2)              % WURZEL-COSINUS-ROLL-OFF oder
  %                                         COSINUS-ROLL-OFF
  rolloff = imp_par;

  if (rolloff < 0) | (rolloff > 1)
    error(['ERROR (get_impulse.m):  Roll-off factor must be in the ',...
           'range 0<=rolloff<=1!']);
  end;

  if (imp_no==1)                          % #####  S2a:  WURZEL-COSINUS-
    %                                                          ROLL-OFF-IMPULS
    roll4 = 4*rolloff;                    % ===== Hilfsgroessen berechnen
    x0 = roll4*tn_pos;                    % x0==0  fuer rolloff=0
    x1 = pi*(1+rolloff)*tn_pos;           % x1~=0, da   tn_pos > 0
    x2 = pi*(1-rolloff)*tn_pos;           % x2==0  fuer rolloff=1
    x3 = 1-x0.*x0;                        % x3==0  fuer tn_pos = 1/roll4
    x4 = pi*tn_pos/2;                     % x4~=0, da   tn_pos > 0

    if roll4 ~= 0                         % ===== "Divide by zero" warning
      %                                     vermeiden
      index = find( (abs(tn_pos-1/roll4) < eps) );
      if length(index) > 0                % tn_pos = 1/roll4  =>  x3=0
        x3(index) = eps;                  % set the corresponding value of x3
        %                                   to about 2e-16
      end;
    end;                                  % ===== WURZEL-COSINUS-ROLL-OFF FORMEL
    %                                       fuer pos. Zeiten
    %                                       fN: Nyquistfrequenz.  Wert so 
    %                                           manipuliert, dass
    fN = 1/(2*(1+rolloff*(4/pi-1)));      %     der Hauptkoeffizient gleich 
    %                                           eins ist (= fN/(A8))
    gn_pos = fN*(x0.*cos(x1)+sin(x2)) ./ (x3.*x4);  % s. Gleichung (A10)

    if roll4 ~= 0                         % ===== Werte an kritischen pos.
      %                                     Zeiten richtig setzen
      if length(index) > 0                % Sonderbehandlung von
        %                                   tn_pos = 1/roll4, Gl. (A9)
        gn_pos(index) = -2*fN*rolloff*( 2/pi*cos(pi*(1+rolloff)/roll4) ...
                                       - cos(pi*(1-rolloff)/roll4) );
      end;
    end;
    gn = [flipud(gn_pos); 1; gn_pos];     % Hauptkoeffizient liegt jetzt bei
    %                                       LT*M_over+1!
    if norm_f==1, gn=gn/(2*fN*M_over); end;% Normierung, so dass Max{|Gn|}=1,
    %                                        wobei Gn:=FFT{gn}
    imp_str = sprintf('SqrtCosRoll (r=%3.1f)', rolloff);

  else                                    % #####  S2b:  COSINUS-ROLL-OFF-IMPULS
    roll2 = 2*rolloff;                    % ===== Hilfsgroessen berechnen
    x0 = tn_pos;                          % x1~=0, da   tn_pos > 0
    x1 = rolloff*pi*tn_pos;               % x1==0  fuer rolloff=0
    x2 = roll2*tn_pos;
    x2 = 1-x2.*x2;                        % x2==0  fuer tn_pos = 1/roll2

    if roll2 ~= 0                         % ===== "Divide by zero" warning
      %                                     vermeiden
      index = find( (abs(tn_pos-1/roll2) < eps) );
      if length(index) > 0                % tn_pos = 1/roll2  =>  x2=0
        x2(index) = eps;                  % set the corresponding value of
        %                                   x2 to about 2e-16
      end;
    end;                                  % ===== COSINUS-ROLL-OFF FORMEL 
    %                                       fuer pos. Zeiten
    %                                       fN: Nyquistfrequenz.  Wert so
    %                                           manipuliert, dass 
    fN = 1/2;                             %     der Hauptkoeffizient gleich
    %                                           eins ist (= fN/(A8))
    gn_pos = 2*fN*sinc(x0).*cos(x1)./x2;  % s. Gleichung (2.1.20b) mit 
    %                                       \omega_N t = \pi tn_pos

    if roll2 ~= 0                         % ===== Werte an kritischen pos.
      %                                     Zeiten richtig setzen
      if length(index) > 0                % Sonderbehandlung von
        %                                   tn_pos = 1/roll2, Gl. (A9)
        gn_pos(index) = 2*fN*sinc(1/roll2)*pi/4;
      end;
    end;
    gn = [flipud(gn_pos); 1; gn_pos];     % Hauptkoeffizient liegt jetzt bei
    %                                       LT*M_over+1!
    if norm_f==1, gn=gn/M_over; end;      % Normierung, so dass Max{|Gn|}=1,
    %                                       wobei Gn:=FFT{gn}
    imp_str = sprintf('CosRoll (r=%3.1f)', rolloff);
  end;


elseif (imp_no==3) | (imp_no==4)          % GMSK oder GAUSS
  f3dB_T = imp_par;

  if f3dB_T < 0
    error('ERROR (gauss.m):  f3dB_T must be positive!');
  end;

  if (imp_no==3)                          % #####  S2c:  GMSK-IMPULS
    if f3dB_T==0
      gn_pos = ones(size(tn_pos));        % ===== "Divide by zero" warning
      %                                     vermeiden
      gn = [flipud(gn_pos); 1; gn_pos];   % Hauptkoeffizient liegt jetzt bei
      %                                     LT*M_over+1!
      if norm_f==1
        gn = gn/length(gn);               % Normierung, so dass Max{|Gn|}=1,
        %                                   wobei Gn:=FFT{gn}
      end;
    else
      alpha = sqrt(2/log(2))*pi*f3dB_T;   % ===== GAUSS*RECT(t/T)-FORMEL
      %                                     fuer pos. Zeiten
      gn_pos = 0.5*(erf(alpha*(tn_pos+0.5))-erf(alpha*(tn_pos-0.5)));
      gn_pos = gn_pos/erf(alpha/2);       % Normierung der Amplitude auf 1
      gn = [flipud(gn_pos); 1; gn_pos];   % Hauptkoeffizient liegt jetzt bei
      %                                     LT*M_over+1!
      if norm_f==1
        gn = gn*erf(alpha/2)/M_over;      % Normierung, so dass Max{|Gn|}=1,
        %                                   wobei Gn:=FFT{gn}
      end;
    end;
    imp_str = sprintf('GMSK (f3T=%4.2f)', f3dB_T);

  else                                    % #####  S2d:  GAUSS-IMPULS
    gn_pos = exp(-2/log(2)*(pi*f3dB_T*tn_pos).^2);  % ===== Pure GAUSS-FORMEL
    %                                                 fuer pos. Zeiten
    gn = [flipud(gn_pos); 1; gn_pos];     % Hauptkoeffizient liegt jetzt bei
    %                                       LT*M_over+1!
    if norm_f==1                          % Normierung, so dass Max{|Gn|}=1,
      %                                     wobei Gn:=FFT{gn}
      gn = gn*sqrt(2*pi/log(2))*f3dB_T/M_over;
    end;
    imp_str = sprintf('Gauss (f3T=%4.2f)', f3dB_T);
  end;
end;
% #####  S3:  Impulsantwort fuer t<0 gerade ergaenzen
if nargout>1
  tn = [-flipud(tn_pos); 0; tn_pos];
end;
if nargout>2
  imp_name = imp_str;
end;

% ### EOF ######################################################################
