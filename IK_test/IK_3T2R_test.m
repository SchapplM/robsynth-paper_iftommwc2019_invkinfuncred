% Teste die inverse Kinematik mit 3T2R mit verschiedenen Systemen

clear
clc

RobotNames = {'kuka6dof', 'kuka5dof', 'S_UPS1', 'lwr4p'};

%% Alle Robotermodelle durchgehen
for mdlname2 = RobotNames
  mdlname = mdlname2{1};
  
  eval(sprintf('TSS = %s_varpar_testfunctions_parameter();', mdlname));
  
  %% Klasse für seriellen Roboter erstellen
  Par_struct = struct('alpha', TSS.alpha, 'a', TSS.a, ...
                      'theta', TSS.theta, 'd', TSS.d, ...
                      'sigma', TSS.sigma, ...
                      'pkin', TSS.pkin, ...
                      'm', TSS.m, 'mrSges', TSS.mrSges, 'Ifges', TSS.Ifges, ...
                      'NJ', TSS.NJ, 'NL', TSS.NL, 'NQJ', TSS.NQJ);
  RS = SerRob(Par_struct, mdlname);
  RS = RS.fill_fcn_handles();
  
  %% Init
  q = TSS.Q(1,:)';
  
  %% Funktionen testweise aufrufen
  
  %% Normale Inverse Kinematik prüfen
  for m = 1:2 % Nach zwei Methoden prüfen
    if RS.NJ >= 6 % Normale IK geht nur mit 6 Roboter-FG oder mehr
      for i = 1:size(TSS.Q,1)
        q = TSS.Q(i,:)';
        T_E = RS.fkineEE(q);
        xE = [T_E(1:3,4); r2rpy(T_E(1:3,1:3))];
        q0 = q-10*pi/180*(0.5-rand(RS.NQJ,1)); % Anfangswinkel 20° neben der Endstellung
        if m == 1
          q_test = RS.invkin1(xE, q0);
        else
          q_test = RS.invkin2(xE, q0);
        end
        T_E_test = RS.fkineEE(q_test);
        test_T = T_E\T_E_test - eye(4);
        % test_q = q-q_test
        % [q, q_test]
        if any(abs(test_T(:)) > 1e-10)
          error('DK/IK stimmt nicht');
        end
      end
      fprintf('%s: Inverse Kinematik Variante %d getestet\n', mdlname, m);
    else
      fprintf('%s: Inverse Kinematik Variante %d nicht getestet\n', mdlname, m);
    end
  end
  %% Prüfe IK mit Aufgabenredundanz
  n_iO = 0;
  for i = 1:size(TSS.Q,1)
    
    q = TSS.Q(i,:)';
    T_E = RS.fkineEE(q);
    xE = [T_E(1:3,4); r2rpy(T_E(1:3,1:3))];
    xE(6) = 0; % Rotation um z-Achse des EE interessiert nicht.
    q0 = q-20*pi/180*(0.5-rand(RS.NQJ,1)); % Anfangswinkel 20° neben der Endstellung
    T_E0 = RS.fkineEE(q0);
    q_test = RS.invkin2(xE, q0, true);
    T_E_test = RS.fkineEE(q_test);
    test_T = T_E\T_E_test - eye(4);
    test_T = test_T(:,[3,4]); % Spalten mit x-y-Einheitsvektoren lassen sich nicht vergleichen.
    if any(abs(test_T(:)) > 1e-10)
      % Teilweise konvergiert die IK nicht, wenn der Abstand zu groß ist.
      warning('DK/IK stimmt nicht für Aufgabenredundanz');
    else
      n_iO = n_iO+1;
    end
  end
  fprintf('%s: Inverse Kinematik Variante 2 mit Aufgabenredundanz=1 getestet. %d/%d erfolgreich\n', ...
    mdlname, n_iO, size(TSS.Q,1));
  
  %% Vergleich Zwangsbedingungen mit selbst berechneten Winkeln
  % Diese Berechnung ist Roboterunabhängig
  for i = 2:size(TSS.Q,1)
    % Rahmenbedingungen für folgende Berechnungen:
    % valide direkte Kinematik berechnen für Soll-Pose xE
    qx = TSS.Q(i-1,:)';
    T_Ex = RS.fkineEE(qx);
    xE = [T_Ex(1:3,4); r2rpy(T_Ex(1:3,1:3))];
    
    % zufälliger anderer Winkel für Aufstellung der Zwangsbedingungen
    q = TSS.Q(i,:)';
    T_E = RS.fkineEE(q);
    R_0_Eq = T_E(1:3,1:3);
    
    % Zwangsbedingungen auf Weg 1 (volle Euler-Winkelkonvention)
    Phi1 = RS.constr2(q, xE);
    %% Versuch 1: Nachvollziehen der vollen Transformation
    % 0 -> TA -> Ex -> Eq -> 0
    % ZB auf Weg 2 (nur xy-Winkel auf beiden Wegen)
    % Bild 3, Weg 0 -> TA (über beta1,beta2)
    R_0_TA = rotx(xE(4)) * roty(xE(5));
    % Bild 3, Weg Ex -> TA -> 0 -> Eq (über beta3,beta2,beta1)
    R_Ex_Eq = (R_0_TA*rotz(xE(6)))' * R_0_Eq;
    % manuelle Berechnung der Winkel der ZYX-Euler-Notation
    % Siehe Aufzeichnungen vom 14.08.2018
    alpha1=atan2(R_Ex_Eq(3,2), R_Ex_Eq(3,3)); % 14.8., Gl. 20
    alpha2=atan2(-R_Ex_Eq(3,1), sqrt(R_Ex_Eq(1,1)^2+R_Ex_Eq(2,1)^2)); % 14.8., Gl. 22
    alpha3=atan2(R_Ex_Eq(2,1),R_Ex_Eq(1,1)); % 14.8., Gl. 23
    % ZB selbst nachgerechnet (Weg 2)
    Phi2 = [Phi1(1:3); alpha1; alpha2; alpha3];
    % Bild 3, Weg 0 -> Eq -> Ex (über alpha1,alpha2,alpha3+beta3)
    R_0_TA_test = R_0_Eq * (rotz(alpha3+xE(6)) * roty(alpha2) * rotx(alpha1))';
    
    % Vergleich der Rotationsmatrix auf zwei Wegen berechnet
    R_test = R_0_TA - R_0_TA_test;
    if any(abs(R_test(:)) > 1e-10)
      error('Die Rotationsmatrix stimmt nicht');
    end
    Phi_test = Phi1 - Phi2;
    if any(abs(Phi_test) > 1e-10)
      error('Die Winkel alpha stimmen nicht');
    end
%     continue
    %% Versuch 2: Nachvollziehen der Transformation nur über Symm.-Achse
    % TODO: Das funktioniert noch nicht!
    % 0 -> TA -> Eq -> 0
    % ZB auf Weg 2 (nur xy-Winkel auf beiden Wegen)
    % Bild 3, Weg 0 -> TA (über beta1,beta2)
    R_0_TA = rotx(xE(4)) * roty(xE(5));
    % Bild 3, Weg TA -> Ex (über beta3)
    R_0_Ex = R_0_TA*rotz(xE(6));
    % Bild 3, Weg TA -> 0 -> Eq (über beta1,beta2)
    R_Ex_Eq = (R_0_Ex)' * R_0_Eq; % *rotz(xE(6)
    
    % manuelle Berechnung der Winkel der YX-Euler-Notation
    % alpha1,alpha2,alpha3 siehe vorheriger Abschnitt (aus R_Ex_Eq)
    % ZB selbst nachgerechnet (Weg 2)
    Phi2 = [Phi1(1:3); alpha1; alpha2; NaN];
    
    % Bild 3, Weg 0 -> Eq -> TA (über alpha2,alpha1)
    R_0_TA_test = R_0_Eq * (roty(alpha2) * rotx(alpha1))';
    R_0_Eq_test = R_0_TA * (roty(alpha2) * rotx(alpha1));
    % R_0_Eq_test - R_0_Eq;
    R_test = R_0_TA - R_0_TA_test;
    if any(abs(R_test(:,3)) > 1e-10)
      error('Die z-Achse stimmt nicht');
    end
    Phi_test = Phi1 - Phi2;
    if any(abs(Phi_test) > 1e-10)
      error('Die Winkel alpha stimmen nicht');
    end
    delta_phiz = r2rpy(R_0_TA \ R_0_TA_test);
%     normalize_angle(+xE(6)+Phi1(6))
  end
  %% TODO: Prüfe IK mit Aufgabenredundanz und Nebenbedingungen
  
  %% TODO: Prüfe IK mit Aufgabenredundanz für Trajektorie
end