      Subroutine Move(Av1,Av2,Delta)
      Implicit None

      Include 'commons.inc'

C     Displace A Randomly Selected Particle

      Logical Laccept
      Integer Ib,Idi,Ipart,Selectint
      Double Precision Rxn,Ryn,Rzn,Xi,Yi,Zi,Ran_Uniform,Unew,Uold
     $     ,Virnew,Virold,Av1,Av2,Delta

      Ipart = Selectint(Npart)
      Ib    = Ibox(Ipart)
      Idi   = Id(Ipart)

      Rxn = Rx(Ipart) + (2.0d0*Ran_Uniform()-1.0d0)*Delta
      Ryn = Ry(Ipart) + (2.0d0*Ran_Uniform()-1.0d0)*Delta
      Rzn = Rz(Ipart) + (2.0d0*Ran_Uniform()-1.0d0)*Delta

C     Put Back In The Box

      If(Rxn.Lt.0.0d0) Then
         Rxn = Rxn + Box(Ib)
      Elseif(Rxn.Gt.Box(Ib)) Then
         Rxn = Rxn - Box(Ib)
      Endif

      If(Ryn.Lt.0.0d0) Then
         Ryn = Ryn + Box(Ib)
      Elseif(Ryn.Gt.Box(Ib)) Then
         Ryn = Ryn - Box(Ib)
      Endif

      If(Rzn.Lt.0.0d0) Then
         Rzn = Rzn + Box(Ib)
      Elseif(Rzn.Gt.Box(Ib)) Then
         Rzn = Rzn - Box(Ib)
      Endif

      Xi = Rx(Ipart)
      Yi = Ry(Ipart)
      Zi = Rz(Ipart)

      Call Epart(Ib,Virold,Uold,Xi,Yi,Zi,Ipart,Idi)
      Call Epart(Ib,Virnew,Unew,Rxn,Ryn,Rzn,Ipart,Idi)

      Call Accept(Dexp(-Beta*(Unew-Uold)),Laccept)

      Av2 = Av2 + 1.0d0

C     Accept Or Reject

      If(Laccept) Then
         Av1 = Av1 + 1.0d0

         Etotal(Ib) = Etotal(Ib) + Unew   - Uold
         Vtotal(Ib) = Vtotal(Ib) + Virnew - Virold

         Rx(Ipart) = Rxn
         Ry(Ipart) = Ryn
         Rz(Ipart) = Rzn
      Endif

      Return
      End
