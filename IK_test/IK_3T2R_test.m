% Teste die inverse Kinematik mit 3T2R mit verschiedenen Systemen

clear
clc

RobotNames = {'kuka5dof', 'S_UPS1'};

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
  
  %% Inverse Kinematik prüfen
  % Normale Inverse Kinematik

  for m = 1:2 % Nach zwei Methoden prüfen
    if RS.NJ >= 6 % Normale IK geht nur mit 6 Roboter-FG oder mehr
      for i = 1:size(TSS.Q,1)
        q = TSS.Q(i,:)';
        T_E = RS.fkineEE(q);
        xE = [T_E(1:3,4); r2rpy(T_E(1:3,1:3))];
        q0 = rand(RS.NQJ,1);
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
  n_iO = 0;
  for i = 1:size(TSS.Q,1)
    % Prüfe mit Aufgabenredundanz
    q = TSS.Q(i,:)';
    T_E = RS.fkineEE(q);
    xE = [T_E(1:3,4); r2rpy(T_E(1:3,1:3))];
    xE(6) = 0; % Rotation um z-Achse des EE interessiert nicht.
    q0 = q-0.1*(0.5-rand(RS.NQJ,1));
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
end