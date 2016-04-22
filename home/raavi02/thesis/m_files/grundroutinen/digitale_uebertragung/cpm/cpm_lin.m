% ##############################################################################
% ##  cpm_lin.m : Berechnung der Elementarimpulse cµ(iT)                      ##
% ##############################################################################
%
% Aufruf:   [ci, tTci, qCPMi, alpha_muenue] = cpm_lin(max_mu, struct_cpm_ini);
%
% Eingabe:  max_mu = Anzahl der zu beruecksichtigenden Elementarimpulse c_mu(i)
%                    [defautl: max_mu = 0]
%           struct_cpm_ini = CPM-Parameter aus cpm_ini
%
% Ausgabe:  ci           = Grund-Elementarimpulse im erhoehten Symboltakt iT/w
%                          der Laenge (2*LT+1)*w+1 fuer gerades und
%                          (2*LT+1)*w   fuer ungerades w
%                          N.B : max_mu = 0: ci ist ein Spaltenvektor
%                          max_mu > 0: ci ist eine Matrix, in der alle
%                          berechneten Grund-Elementarimpulse in den jeweiligen
%                          Spalten aufgefuehrt sind
%           tTci         = auf die Symboldauer T normierte Zeitachse im
%                          erhoehten Symtakt als Spaltenvektor bezogen auf den
%                          Grund-Elementarimpuls ci
%           qCPMi        = integrierter Elementarimpuls
%           alpha_muenue = Matrix bzw. Vektor (abhaengig von max_mu) des
%                          binaerer Indizes
%
% Benoetigte Routinen:  cpm_ini, get_impulse
%
% Quelle:    - [Kam92] K.D. Kammeyer, Nachrichtenuebertragung, B.G.
%                      Teubner Stuttgart 1992, Kapitel 11
%            - [Ben95] Marcus Benthin, Technische Berichte, Blaupunkt,
%                      Okt 1995 Nr. 1
%
% Thorsten Petermann                           31.03.2000 (last update)

function [ci, tTci, qCPMi, alpha_muenue] = cpm_lin(max_mu, struct_cpm_ini);

if nargin<1               % Setzen/Auslesen der Simulationsparameter
  struct_cpm_ini = struct('L', 4, 'impulsform', 'gmsk', 'eta', 0.5, 'M', M,...
                          'w', 4, 'f3dBT', 0.3); max_mu = 0;
else
  if strcmp('GAUSS',struct_cpm_ini.impulsform)
    struct_cpm_ini.impulsform = 'gmsk';
  end
end

L = struct_cpm_ini.L;
impulsform = struct_cpm_ini.impulsform;
eta = struct_cpm_ini.eta;
M = struct_cpm_ini.M;
w = struct_cpm_ini.w;
f3dBT = struct_cpm_ini.f3dBT;

LT    = ceil(L/2);                     % einseitige Beobachtungslaenge
gCPM  = get_impulse (impulsform, f3dBT, 512, LT, 1, 1);
% 512-fach ueberabgetasteter Elementarimpuls
qCPM  = cumsum(gCPM);                  % Integration
qCPMi = qCPM(1:512/w:length(qCPM));    % Downsampling auf w-Takt
clear qCPM;                            % Loeschen des 512-fach ueberabgetast.
%                                        Elementarimpulses

L_qCPMi = length(qCPMi);               % Laenge des CPM-Elementarimpulses
%                                        im w-Takt
LT = floor((L_qCPMi-1)/(2*w));         % Berechnung der einseitigen
%                                        Beobachtungslaenge

if max_mu > 2^(2*LT-1)-1               % sollte die Anzahl der gewaehlten
  max_mu = 2^(2*LT-1)-1;               % Approximationen den Wert 2^(2*LT-1)-1
end                                    % ueberschreiten, so  wird diese auf den
%                                        maximal moeglichen Wert gesetzt

psiCPMi = [qCPMi; 1-qCPMi(2:L_qCPMi)]; % Berechnung des Phasenimpulses
%                                        (siehe [Kam92] 12.1.10)
p0i = sin(pi*eta*psiCPMi)/sin(pi*eta); % Berechnung des Grundimpulses
%                                        (siehe [Kam92] 12.1.11a)
L_p0i = length(p0i);                   % Laenge des Grundimpulses

if max_mu == 0                        % Berechnung des Grund-Elementarimpulses
  %                                     c0(t)(siehe [Kam92] 12.1.13)
  ci_tmp = 1;                         % Anfangswert fuer das Grundimpuls-Produkt

  for ii = 0:2*LT-1 % Schleife zur Berechnung der zeitverschobenen Versionen...
    % des Grundimpulses (siehe [Kam92] 12.1.11b)
    pi_tmp = p0i(ii*w+1:L_p0i-(2*LT-1-ii)*w);
    ci_tmp = ci_tmp.*pi_tmp;           % Produkt saemtlicher zeitverschobener
  end                                  % Versionen des Grundimpulses

  ci = ci_tmp(:);                      % Erzeugung eines Spaltenvektors
  ci_left = flipud(ci((2*LT+1)*w+1:-1:1));
  ci_right = ci((2*LT+1)*w+1:length(ci));
  % Berechnung der linken und rechten Haelfte des abgetas-...
  % teten Grund-Elementarimpulses c0(t) vom Maximalwert aus
  ci = [ci_left; ci_right(2:length(ci_right))];
  alpha_muenue = zeros(1,2*LT-1); % Nullvektor als binaerer
  % Index (Zeilenvektor) (siehe [Kam92] 12.1.12)

elseif max_mu >= 1 % Berechnung saemtlicher Elementarimpulse
  % cµ(t) (siehe [Kam92] 12.1.13)
  pi_mat = [];
  ci = [];

  for ii = 0:4*LT-1 % Schleife zur Berechnung der zeitverschobenen
    % Versionen des Grundimpulses (siehe [Kam92] 12.1.11b)
    pi_tmp = [p0i(ii*w+1:L_p0i).' zeros(1,ii*w)];
    pi_tmp = pi_tmp(1:(2*LT+1)*w+1);
    % Reduzierung der Laenge des temporaeren Grundimpulses...
    % auf den relevanten Bereich 2*LT+1 (siehe [Kam92] 12.1.18)
    pi_mat = [pi_mat; pi_tmp];
    % Matrix, die saemtlichen zeitverschobenen Versionen des...
  end                                  % Grundimpulses zeilenweise enthaelt

  h_vec = (1:2^(2*LT-1)).'; % Hilfsvektor zur Umwandlung des dezimalen
  % Indizes µ in einen binaeren Index alpha (siehe [Kam92] 12.1.12)

  for ii = 1:2*LT-1 % Schleife zur Berechnung des binaeren Indizes alpha als...
    % (mue*nue)-Matrix: z.B.:     nue ->
    %                         0 = 0 0 0  mue
    %                         1 = 0 0 1   |
    %                            usw.     v
    alpha_muenue(h_vec,ii) = rem(ceil(h_vec/(2^(2*LT-1-ii))),2);
  end                                  % (siehe [Kam92] 12.1.12)

  alpha_muenue = flipud(fliplr(alpha_muenue));
  % 180°-Drehung der Matrix
  nue = 2:2*LT;                        % Laufvariable: (1 <= nue <= 2*LT-1)+1

  for mue = 0:max_mu  % Schleife ueber alle moeglichen Approximationen...
    % (max: 2^(2*LT-1)) (siehe [Kam92] 12.1.16)
    ind = [1 nue+2*LT*alpha_muenue(mue+1,:)];
    % Erzeugung eines Hilfsindizes zur Ermittlung der rele-...
    % vanten Grundimpulse
    ci_tmp = (prod(pi_mat(ind,:))).';  % Multiplikation der relevanten Zeilen
    % der Grundimpulsmatrix zur Erzeugung der temporaeren Elementarimpulse
    ci_left_tmp = flipud(ci_tmp((2*LT+1)*w+1:-1:1));
    ci_right_tmp = ci_tmp((2*LT+1)*w+1:length(ci_tmp));
    % Berechnung der linken und rechten Haelfte der abgetas-...
    % teten Grund-Elementarimpulse cµ(t) vom Maximalwert aus
    ci_tmp = [ci_left_tmp; ci_right_tmp(2:length(ci_right_tmp))];
    ci = [ci ci_tmp]; % Matrix mit allen berechneten Elementarimpulsen
    %                   (spaltenweise)
  end
end

tTci_left = (fliplr(LT+1/2:-1/w:0)).';
tTci_right = (LT+1/2:1/w:2*LT+1).';
tTci = [tTci_left; tTci_right(2:length(tTci_right))];
% Spaltenvektor der normierten Zeitachse fuer den Grund-...
% elementarimuls ck

% ### EOF ######################################################################
