% Teste die inverse Kinematik mit 3T2R mit verschiedenen Systemen
% Untersuche insbesondere, wie die Eigenschaften der reziproken
% Euler-Winkel sich verhalten

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

clear
clc

%% Benutzereingaben
% Prüfung repräsentativer Roboter

Robots = {{'S4RRPR1', 'S4RRPR1_KUKA1'}, ... % 4DoF SCARA
          {'S5RRRRR1', 'S5RRRRR1_KUKA1'}, ... % 5DoF Paletizer
          {'S6RRRRRR10V2', 'S6RRRRRR10V2_KUKA1'}, ... % 6DoF Industrial Robot
          {'S7RRRRRRR1', 'S7RRRRRRR1_LWR4P'}}; % 7DoF Lightweight Robot
% Folgende Zeilen zur Prüfung einzelner Roboter einkommentieren:
% Robots = {{'S4RRPR1', 'S4RRPR1_KUKA1'}};
% Robots = {{'S5RRRRR1', 'S5RRRRR1_KUKA1'}};
% Robots = {{'S6RRRRRR10V2', 'S6RRRRRR10V2_KUKA1'}};
% Robots = {{'S7RRRRRRR1', 'S7RRRRRRR1_LWR4P'}};

% Einstellungen
use_mex_functions = true; % mit mex geht es etwas schneller, dafür ist debuggen schwieriger
minimal_test = true; % Nur sehr wenige zufällige Winkel testen (geht schneller)
% Endeffektor-Transformation ungleich Null festlegen, um zu prüfen, ob die
% Implementierung davon richtig ist
r_W_E = [0.1;0.2;0.3];
phi_W_E = [20; 40; 50]*pi/180;
%% Alle Robotermodelle durchgehen
for Robot_Data = Robots
  SName = Robot_Data{1}{1};
  RName = Robot_Data{1}{2}; 
  %% Klasse für seriellen Roboter erstellen
  % Instanz der Roboterklasse erstellen
  RS = serroblib_create_robot_class(SName, RName);
  RS.fill_fcn_handles(use_mex_functions, true);
  
  % Grenzen festlegen (für Zusatz-Optimierung)
  RS.qlim = repmat([-pi, pi], RS.NQJ, 1);
  RS.update_EE(r_W_E, phi_W_E, []);
  
  % Test-Einstellungen generieren
  TSS = RS.gen_testsettings(true, false);
  if minimal_test
    TSS.Q = TSS.Q(1:3,:);
  end
  fprintf('Starte Untersuchung für %s\n', RS.descr);
  
  %% Init 
  % Bereich reduzieren: Bei Addition von Zufallswinkeln darf nicht pi
  % überschritten werden.
  TSS.Q(abs(TSS.Q(:))>150*pi/180) = 0;
  
  %% Funktionen testweise aufrufen
  q0 = TSS.Q(1,:)';
  RS.fkine(q0);
  RS.jacobig(q0);
  %% Normale Inverse Kinematik prüfen
  for m = 1:2 % Nach zwei Methoden prüfen
    if RS.NJ >= 6 % Normale IK geht nur mit 6 Roboter-FG oder mehr
      for i = 1:size(TSS.Q,1)
        q = TSS.Q(i,:)';
        T_E = RS.fkineEE(q);
        xE = [T_E(1:3,4); r2eulxyz(T_E(1:3,1:3))];
        q0 = q-10*pi/180*(0.5-rand(RS.NQJ,1)); % Anfangswinkel 20° neben der Endstellung
        q_test = RS.invkin(xE, q0, struct('constr_m', m));
        T_E_test = RS.fkineEE(q_test);
        test_T = T_E\T_E_test - eye(4);
        % test_q = q-q_test
        % [q, q_test]
        if any(abs(test_T(:)) > 1e-9)
          error('DK/IK stimmt nicht');
        end
      end
      fprintf('%s: Inverse Kinematik Variante %d getestet\n', SName, m);
    else
      fprintf('%s: Inverse Kinematik Variante %d nicht getestet\n', SName, m);
    end
  end
  
  %% Aufgabenredundanz: Prüfe, ob die Jacobi-Matrix durch beta3 beeinflusst wird
  % Entspricht Prüfung von Gl. 18 und Gl. 22
  for i = 1:size(TSS.Q,1)
    q = TSS.Q(i,:)';
    xE_soll = rand(6,1);
    for j = 1:10
      xE_soll(6) = rand(); % zufälliger neuer EE-Winkel
      % Prüfe Jacobi-Matrix
      dpq=RS.constr2grad_rq(q, xE_soll, true);
      dqp_IK = dpq(1:2,:);
      if j > 2
        if any(abs(dqp_IK_alt(:) - dqp_IK(:)) > 1e-10)
          error('Jacobi-Matrix für Aufgabenredundanz hat sich durch anderen Winkel beta3 geändert. Darf nicht sein.');
        end
      end
      dqp_IK_alt = dqp_IK;
      % Prüfe ZB
      Phi = RS.constr2(q, xE_soll, true);
      Phi_IK = Phi(5:6,:);
      if j > 2
        if any(abs(Phi_IK(:) - Phi_IK_alt(:)) > 1e-10)
          error('Rotatorische Zwangsbedingungen für Aufgabenredundanz hat sich durch anderen Winkel beta3 geändert. Darf nicht sein.');
        end
      end
      Phi_IK_alt = Phi_IK;
    end
  end
  %% Prüfe IK mit Aufgabenredundanz
  % Siehe Sec. 3.1 im Paper
  n_iO = 0;
  for i = 1:size(TSS.Q,1)
    warning off
    q = TSS.Q(i,:)';
    T_E = RS.fkineEE(q);
    xE = [T_E(1:3,4); r2eulxyz(T_E(1:3,1:3))];
    xE(6) = 0; % Rotation um z-Achse des EE interessiert nicht.
    q0 = q-20*pi/180*(0.5-rand(RS.NQJ,1)); % Anfangswinkel 20° neben der Endstellung
    T_E0 = RS.fkineEE(q0);
    q_test = RS.invkin(xE, q0, struct('constr_m', 2, 'task_red', true));
    T_E_test = RS.fkineEE(q_test);
    test_T = T_E\T_E_test - eye(4);
    test_T = test_T(:,[3,4]); % Spalten mit x-y-Einheitsvektoren lassen sich nicht vergleichen.
    if any(abs(test_T(:)) > 1e-9) || any(isnan(test_T(:)))
      % Teilweise konvergiert die IK nicht, wenn der Abstand zu groß ist.
      warning on
      warning('DK/IK stimmt nicht für Aufgabenredundanz. Delta_x = %1.5e, Delta_z = %1.5e', norm(test_T(:,2)), norm(test_T(:,1)));
    else
      n_iO = n_iO+1;
    end
  end
  fprintf('%s: Inverse Kinematik Variante 2 mit Aufgabenredundanz=1 getestet. %d/%d erfolgreich\n', ...
    SName, n_iO, size(TSS.Q,1));
  
  %% Vergleich Zwangsbedingungen mit selbst berechneten Winkeln
  % Diese Berechnung ist Roboterunabhängig
  for i = 2:size(TSS.Q,1)
    % Rahmenbedingungen für folgende Berechnungen:
    % valide direkte Kinematik berechnen für Soll-Pose xE
    qx = TSS.Q(i-1,:)';
    T_Ex = RS.fkineEE(qx);
    xE = [T_Ex(1:3,4); r2eulxyz(T_Ex(1:3,1:3))];
    
    % zufälliger anderer Winkel für Aufstellung der Zwangsbedingungen
    q = TSS.Q(i,:)';
    T_E = RS.fkineEE(q);
    R_0_Eq = T_E(1:3,1:3);
    
    % Zwangsbedingungen auf Weg 1 (volle Euler-Winkelkonvention)
    Phi1 = RS.constr2(q, xE, true);
    %% Versuch 1: Nachvollziehen der vollen Transformation
    % Siehe Fig. 1 im Paper. Berechnung des Residuums Phi mit der Formel
    % 0 -> TA -> Ex -> Eq -> 0
    % ZB auf Weg 2 (nur xy-Winkel auf beiden Wegen)
    % Bild 3, Weg 0 -> TA (über beta1,beta2)
    R_0_TA = rotx(xE(4)) * roty(xE(5));
    % Bild 3, Weg Ex -> TA -> 0 -> Eq (über beta3,beta2,beta1)
    R_Ex_Eq = (R_0_TA*rotz(xE(6)))' * R_0_Eq;
    % manuelle Berechnung der Winkel der ZYX-Euler-Notation
    % Siehe Aufzeichnungen vom 14.08.2018 (dort extrinsische Drehungen, hier intrisisch)
    % (daher hier andere Reihenfolge der alpha)
    alpha1=atan2(R_Ex_Eq(3,2), R_Ex_Eq(3,3)); % 14.8., Gl. 20
    alpha2=atan2(-R_Ex_Eq(3,1), sqrt(R_Ex_Eq(1,1)^2+R_Ex_Eq(2,1)^2)); % 14.8., Gl. 22
    alpha3=atan2(R_Ex_Eq(2,1),R_Ex_Eq(1,1)); % 14.8., Gl. 23
    % ZB selbst nachgerechnet (Weg 2)
    Phi2 = [Phi1(1:3); alpha3; alpha2; alpha1];
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

    %% Versuch 2: Nachvollziehen der Transformation nur über Symm.-Achse
    % Siehe Fig. 1 im Paper. Händisches berechnen der Transformation aus
    % dem Bild.
    % ZB auf Weg 2 (nur xy-Winkel auf alpha-Weg)
    % Paper Fig. 1, Weg 0 -> TA1 (über beta1,beta2)
    R_0_TA1 = rotx(xE(4)) * roty(xE(5));
    
    R_TA1_Eq = (R_0_TA1*rotz(rand))' * R_0_Eq;
    % manuelle Berechnung der Winkel der ZYX-Euler-Notation
    % Siehe Aufzeichnungen vom 14.08.2018
    % Die Winkel alpha sind die gleichen, werden aber anders berechnet
    alpha1=atan2(R_TA1_Eq(3,2), R_TA1_Eq(3,3)); % 14.8., Gl. 20
    alpha2=atan2(-R_TA1_Eq(3,1), sqrt(R_TA1_Eq(1,1)^2+R_TA1_Eq(2,1)^2)); % 14.8., Gl. 22

    Phi2 = [Phi1(1:3); NaN; alpha2; alpha1];
    
    % Paper Fig. 1, Weg 0 -> Eq -> TA2 (über alpha2,alpha1)
    R_0_TA2 = R_0_Eq * (roty(alpha2) * rotx(alpha1))';
    R_TA_test = R_0_TA1 - R_0_TA2;
    if any(abs(R_TA_test(:,3)) > 1e-10)
      error('Die z-Achse der Rotationsmatrix TA2 stimmt nicht');
    end

    % Paper Fig. 1, Weg TA -> TA2 (über beta3,alpha3)
    R_TA1_TA2 = rotz(xE(6)+alpha3);
    R_0_TA1_test = R_0_TA2*R_TA1_TA2';
    R_TA1_test = R_0_TA1 - R_0_TA1_test;
    if any(abs(R_TA1_test(:)) > 1e-10)
      error('Die Rotationsmatrix TA stimmt nicht');
    end
    Phi_test = Phi1 - Phi2;
    if any(abs(Phi_test) > 1e-10)
      error('Die Winkel alpha stimmen nicht');
    end
  end
end