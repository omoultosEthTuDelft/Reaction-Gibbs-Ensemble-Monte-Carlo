      Integer Maxpart,Maxcomp,Maxreact
      
      Parameter (Maxpart  = 5000)
      Parameter (Maxcomp  = 2)
      Parameter (Maxreact = 1)

      Double Precision Rgas, to_kPa
      
      Parameter (Rgas     = 8.314d-3)
      Parameter (to_kPa   = 13807d0)

      Integer Nibox,Npart,Ncomp,Nreact,Nstoi(Maxcomp,Maxreact),Ibox(Maxpart)
     $     ,Id(Maxpart),Npbox(2),Nid(Maxcomp,2)

      Double Precision Rx(Maxpart),Ry(Maxpart),Rz(Maxpart),Box(2)
     $     ,lnQ(Maxcomp),Beta,Etotal(2),Vtotal(2),Rcutsq,Rcut2(Maxcomp
     $     ,Maxcomp),Sig2(Maxcomp,Maxcomp),Ecut(Maxcomp,Maxcomp)
     $     ,Eps(Maxcomp,Maxcomp),Press

      Common /Mycomms/ Rx,Ry,Rz,Box,lnQ,Beta,Etotal,Vtotal,Rcutsq,Rcut2
     $     ,Sig2,Ecut,Eps,Press,Nibox,Npart,Ncomp,Nreact,Nstoi,Ibox,Id,Npbox
     $     ,Nid
