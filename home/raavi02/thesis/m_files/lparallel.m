% ##############################################################################
% ##  lparallel.m : Bewerkstelligt die l-fache Parallelschaltung von Codes    ##
% ##############################################################################
%
% Funktionsaufruf: A_out = lparallel(A_in,l)
%------------------------------------------------------
% Eingangsparameter:
%     A_in: IOWEF des Originalcodes
%     l   : Anzahl paralleler Codes
%
% Ausgangsparameter:
%     A_out: IOWEF des l-fachen Codes
% ------------------------------------------------------------------------------
% Volker Kuehn, 20.02.01
% Juergen Rinas, 05.03.01 verwendet jetzt conv2()
% ------------------------------------------------------------------------------

function A_out= lparallel(A_in,l)

A_out=1;
for i=1:l
  A_out=conv2(A_out,A_in);
end;

% ### EOF ######################################################################
