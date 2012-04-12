      FUNCTION TAUISS(d, sm, nu)
c
c calculates the pulse broadening time in ms
c from distance, scattering measure, and radio frequency
c
c input:      d = pulsar distance       (kpc)    
c            sm = scattering measure    (kpc m^{-20/3})
c            nu = radio frequency       (GHz)
c output: tauiss = pulse broadening time (ms) 
c
      real nu
      tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)
      end
c
c
      FUNCTION SCINTBW(d, sm, nu)
c
c calculates the scintillation bandwidth in kHz 
c from distance, scattering measure, and radio frequency
c
c input:        d = pulsar distance       (kpc)    
c              sm = scattering measure    (kpc m^{-20/3})
c              nu = radio frequency       (GHz)
c output: scintbw = scintillation bandwidth (kHz)
c
      real nu
      tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)
      scintbw = 1. / (2. * 3.14159 * tauiss)
      end
c
      FUNCTION SCINTIME(sm, nu, vperp)
c
c calculates the scintillation speed for given distance, galactic
c longitude and latitude, frequency, and transverse velocity      
c
c input:   sm = scattering measure	(kpc m^{-20/3})
c          nu = radio frequency 	(GHz)
c       vperp = psr transverse speed  	(km/s)  
c
c output: scintime = scintillation time (sec)
c 
c usage: should be called with sm = smtau for appropriate
c        line of sight weighting
c reference: eqn (46) of Cordes & Lazio 1991, ApJ, 376, 123.
c
      real nu
      scintime = 2.3 * nu**1.2 * sm**(-0.6) * (100./vperp)
      end
c
c
      FUNCTION THETA_XGAL(sm, nu)
c
c calculates angular broadening for an extragalactic
c source of plane waves
c
c sm = scattering measure
c nu = radio frequency
c theta_xgal = angular broadening FWHM (mas)
c
      real nu
      theta_xgal = 128. * sm**0.6 * nu**(-2.2)
      end
c
      FUNCTION THETA_GAL(sm, nu)
c
c calculates angular broadening for a galactic
c source of spherical waves
c
c sm = scattering measure
c nu = radio frequency
c theta_gal = angular broadening FWHM (mas)
c
      real nu
      theta_gal = 71. * sm**0.6 * nu**(-2.2)
      end
c
      FUNCTION EM (sm)
c
c units of sm are kpc m^{-20/3}
c units of em are pc cm^{-6}      
c
c calculates the emission measure from the scattering measure
c using an assumed outer scale and spectral index of the
c wavenumber spectrum.
c
c for a wavenumber spectrum P_n(q) = q^{-alpha} from q_0 to q_1
c the mean square electron density is
c
c <n_e^2> =~  4pi*[C_n^2 / (alpha - 3) ] * q_0^{3 - alpha)
c
c ( an approximate form that assumes (q_0 / q_1)^{3-alpha} >> 1.
c
c Jim Cordes 18 Dec 1989
c
      data router /1./	! outer scale = 1 pc
      data pc/3.086e+18/
      data alpha/3.6666667/
      data pi/3.14159/
c
      em = sm *
     1         ( (4. * pi * 1000.) / (alpha - 3.) ) * 
     2         (router*pc / (2. * 3.14159) )**(alpha-3.) *
     3         (0.01) ** (20./3.)
c
      return
      end        
