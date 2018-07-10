function xdot = f (x,t)
     
      Shh = 0.058*150 ; # Shh quantity [0,30]
      k_shh = 0.58  ;# dissociation constant shh-ptc bindings [0.58,2.0]
      k_ptc = 8.3*10**-2;  # 1/2maximal concentration of ptc which inhibits smo signlaing
      k_deg = 0.009  ;# degradation constant for all x(1) related proteins
      k_g3rc = 0.012  ;# rate constant for the conversion to signal strenGh
      r_g3b = 1.6*10**-1;  # basal rate of x(2) synthesis
      K_g3rc = 0.1 ; # sensitivity constant of the conversion to signal strenGh
      k_deg_p = 0.09 ; # degradation rate constant for x(4) [0.045,0.071]
      # --------------
      K1 = 8.3*10**-1;
      K2 = 8.3*10**-1;
      c = 1;
      e = 0.5;
      r = 0.2;
      v_max = 2.4*10**-1;
      r_bas = v_max/100;
      v_maxp = 7.5*10**-1;
      r_basp = v_maxp/100;
      xdot(1) = v_max*(((x(2)*K1+x(1)*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(x(2)*K1+x(1)*K2+2*e*x(3)*K1*r)+c**2*(x(2)**2*K1**2+x(1)**2*K2**2+3*e*x(1)*x(3)*K1*K2*r + 3*e**2*x(3)**2*K1**2*r**2 + x(2)*K1*(2*x(1)*K2 + 3*e*x(3)*K1*r))))/(3*c*K1*K2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**2 + c**2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**3 +K1**2*K2**2*(3*x(2)*K1 + 3*x(3)*K1 + (3*x(1) + K1)*K2)))+r_bas*((3*c*K1*K2*(x(2)*K1+ x(1)*K2 + x(3)*K1*r)**2 + c**2*(x(2)*K1 + x(1)*K2 + x(3)*K1*r)**3 + K1**2*K2**2*(3*x(2)*K1 + 3*x(1)*K2 + K1*(K2+ 3*x(3)*r)))/ (3*c*K1*K2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**2 + c**2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**3 + K1**2*K2**2*(3*x(2)*K1 + 3*x(3)*K1 + (3*x(1) +K1)*K2)))-k_deg*x(1);
      xdot(2) = r_g3b/x(4)-x(2)*(k_deg+(k_g3rc/(K_g3rc+((1+Shh/k_shh)/(1+Shh/k_shh+x(4)/k_ptc)))));   
      xdot(3) = x(2)*(k_g3rc/(K_g3rc+(((1+Shh/k_shh)/(1+Shh/k_shh+x(4)/k_ptc)))))-k_deg*x(3)
      xdot(4) = v_maxp*(((x(2)*K1+x(1)*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(x(2)*K1+x(1)*K2+2*e*x(3)*K1*r)+c**2*(x(2)**2*K1**2+x(1)**2*K2**2+3*e*x(1)*x(3)*K1*K2*r + 3*e**2*x(3)**2*K1**2*r**2 + x(2)*K1*(2*x(1)*K2 + 3*e*x(3)*K1*r))))/(3*c*K1*K2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**2 + c**2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**3 +K1**2*K2**2*(3*x(2)*K1 + 3*x(3)*K1 + (3*x(1) + K1)*K2)))+r_basp*((3*c*K1*K2*(x(2)*K1+ x(1)*K2 + x(3)*K1*r)**2 + c**2*(x(2)*K1 + x(1)*K2 + x(3)*K1*r)**3 + K1**2*K2**2*(3*x(2)*K1 + 3*x(1)*K2 + K1*(K2+ 3*x(3)*r)))/ (3*c*K1*K2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**2 + c**2*(x(2)*K1 + x(3)*K1 + x(1)*K2)**3 + K1**2*K2**2*(3*x(2)*K1 + 3*x(3)*K1 + (3*x(1) +K1)*K2)))-k_deg_p*x(4)  
  endfunction

