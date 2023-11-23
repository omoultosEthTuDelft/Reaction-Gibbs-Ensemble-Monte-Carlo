      Subroutine Epart(Ib,Vir,Upot,Xi,Yi,Zi,Ipart,Idp)
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Of Particle Ipart In Box Ib

      Integer I,Ib,Idp,Idi,Ipart
      Double Precision Vir,Upot,Dx,Dy,Dz,R2,Bx,Hbx,Xi,Yi,Zi

      Upot = 0.0d0
      Vir  = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

      Do I=1,Npart

         If(Ibox(I).Eq.Ib.And.I.Ne.Ipart) Then

            Idi = Id(I)
            
            Dx = Rx(I)-Xi
            Dy = Ry(I)-Yi
            Dz = Rz(I)-Zi
            
            If (Dx.Gt.Hbx) Then
               Dx = Dx - Bx
            Elseif (Dx.Lt.-Hbx) Then
               Dx = Dx + Bx
            Endif
                  
            If (Dy.Gt.Hbx) Then
               Dy = Dy - Bx
            Elseif (Dy.Lt.-Hbx) Then
               Dy = Dy + Bx
            Endif
                  
            If (Dz.Gt.Hbx) Then
               Dz = Dz - Bx
            Elseif (Dz.Lt.-Hbx) Then
               Dz = Dz + Bx
            Endif

            R2 = Dx**2 + Dy**2 + Dz**2

            If(R2.Lt.Rcut2(Idi,Idp)) Then
               R2   = Sig2(Idi,Idp)/R2
               R2   = R2*R2*R2
               Upot = Upot + Eps(Idi,Idp)*R2*(R2-1.0d0) - Ecut(Idi,Idp)
               Vir  = Vir  + 12.0d0*Eps(Idi,Idp)*R2*(R2-0.5d0)
            Endif
         Endif
      Enddo

      Return
      End
