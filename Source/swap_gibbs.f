      Subroutine Swap_Gibbs(Av1,Av2)
      Implicit None

C     Particle Swap In The Gibbs Ensemble

      Include 'commons.inc'

      Logical Laccept
      Integer Iadd,Idel,Ipart,Idp,Selectint
      Double Precision Xi,Yi,Zi,Ran_Uniform,Av1,Av2,Upotadd,Upotdel
     $     ,Viradd,Virdel

C     Select At Random Which Box To Add / Remove

      If(Ran_Uniform().Lt.0.5d0) Then
         Iadd = 1
         Idel = 2
      Else
         Iadd = 2
         Idel = 1
      Endif

      Idp = Selectint(Ncomp)
      
      If(Nid(Idp,Idel).Eq.0) Return

C     Select Particle To Be Removed Random From Box Idel

 1    Continue
      Ipart = Selectint(Npart)
      If(Ibox(Ipart).Ne.Idel) Goto 1
      If(Id(Ipart).Ne.Idp) Goto 1

      Xi = Rx(Ipart)
      Yi = Ry(Ipart)
      Zi = Rz(Ipart)

      Call Epart(Idel,Virdel,Upotdel,Xi,Yi,Zi,Ipart,Idp)

      Xi = Ran_Uniform()*Box(Iadd)
      Yi = Ran_Uniform()*Box(Iadd)
      Zi = Ran_Uniform()*Box(Iadd)

      Call Epart(Iadd,Viradd,Upotadd,Xi,Yi,Zi,0,Idp)

      Call Accept((Dble(Nid(Idp,Idel))*(Box(Iadd)**3)/
     &     ((Box(Idel)**3)*Dble(Nid(Idp,Iadd)+1)))*
     &     Dexp(-Beta*(Upotadd-Upotdel)),Laccept)

      Av2 = Av2 + 1.0d0

C     Accept Or Reject

      If(Laccept) Then
         Npbox(Iadd) = Npbox(Iadd) + 1
         Npbox(Idel) = Npbox(Idel) - 1

         Nid(Idp,Iadd) = Nid(Idp,Iadd) + 1
         Nid(Idp,Idel) = Nid(Idp,Idel) - 1
         
         Rx(Ipart) = Xi
         Ry(Ipart) = Yi
         Rz(Ipart) = Zi

         Ibox(Ipart) = Iadd

         Av1 = Av1 + 1.0d0

         Etotal(Iadd) = Etotal(Iadd) + Upotadd
         Etotal(Idel) = Etotal(Idel) - Upotdel

         Vtotal(Iadd) = Vtotal(Iadd) + Viradd
         Vtotal(Idel) = Vtotal(Idel) - Virdel
      Endif
      
      Return
      End
