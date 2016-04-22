% ##############################################################################
% ##  pr_decod.m : Partial-Response-Decodierung                               ##
% ##############################################################################
%
% Aufruf:        b_dach = pr_decod(c,Typ);
%
% Eingabe:       c(i) = Partial-Response-Symbole
%                Typ  = 4 moegliche PR-Codes:
%                         1: Duobinaercode           alpha={1,1}
%                         2: Pseudoternaercode,n=2   alpha={1,-1}
%                         3: Pseudoternaercode,n=3   alpha={1,0,-1}
%                         4: 5-stufiger Code,n=5     alpha={-1,0,2,0,-1}
%
% Ausgabe:       b(i) = decodierte Binaerdaten in 0/1 Schreibweise
%
% Voraussetzung: PR-Coder wurde mit Daten im +/- Format angesteuert
%
% Anmerkung:     keine Unterdrueckung von Ein- und Ausschwingern

function b_dach=pr_decod(c,Typ)

if Typ==1
  b_dach=abs(c)-1;                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  b_dach=-sign(b_dach);              %%%%% Duobinaercode      %%%%%
  b_dach=round((b_dach+1)/2);        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif Typ==2
  b_dach=abs(c)-1;                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  b_dach=sign(b_dach);               %%%%% Pseudoternaer, n=2 %%%%%
  b_dach=round((b_dach+1)/2);        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif Typ==3
  b_dach=abs(c)-1;                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  b_dach=sign(b_dach);               %%%%% Pseudoternaer, n=3 %%%%%
  b_dach=round((b_dach+1)/2);        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif Typ==4
  c=abs(c);                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  b_dach=zeros(1,length(c));         %%%%% 5-stufiger Code    %%%%%
  for i=1:length(c)                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if c(i)>1
      b_dach(i)=1;
    else
    end
    if c(i)>3
      b_dach(i)=0;
    else
    end
  end
else
end

% ### EOF ######################################################################
