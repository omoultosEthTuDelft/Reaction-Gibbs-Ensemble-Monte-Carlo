      Subroutine Energy_Remove(Ib,Uold,Virold,Nremove,Iremove)
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Of Particles That Need To Be Removed
      
      Integer Ib,Nremove,Iremove(Maxpart),I,J,Ii,Jj,Idi,Idj

      Double Precision Uold,Virold,R2,Dx,Dy,Dz,Bx,Hbx

      Logical Lok(Maxpart)

C     Check For The Correct Box
      
      Do I=1,Npart
         If(Ibox(I).Eq.Ib) Then
            Lok(I) = .True.
         Else
            Lok(I) = .False.
         Endif
      Enddo

C     Exclude Interactions Between Multiple Particles That Need To Be Removed

      Do I=1,Nremove
         Lok(Iremove(I)) = .False.
      Enddo

      Uold   = 0.0d0
      Virold = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

C     For Each Particle That Need To Be Removed, Calculate The Interactions With
C     All Other Particles, Except The Ones That Need To Be Removed
      
      Do I=1,Npart

         If(Lok(I)) Then

            Idi = Id(I)

            Do Jj=1,Nremove

               J   = Iremove(Jj)
               Idj = Id(J)
               
               Dx = Rx(I)-Rx(J)
               Dy = Ry(I)-Ry(J)
               Dz = Rz(I)-Rz(J)
            
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

               If(R2.Lt.Rcut2(Idi,Idj)) Then
                  R2   = Sig2(Idi,Idj)/R2
                  R2   = R2*R2*R2
                  Uold = Uold + Eps(Idi,Idj)*R2*(R2-1.0d0) - Ecut(Idi,Idj)
                  
                  Virold = Virold  + 12.0d0*Eps(Idi,Idj)*R2*(R2-0.5d0)
               Endif
            Enddo
         Endif
      Enddo

C     Compute The Energy Between Particles That Need To Be Removed
      
      Do Ii=1,Nremove-1

         I   = Iremove(Ii)
         Idi = Id(I)

         Do Jj=Ii+1,Nremove

            J   = Iremove(Jj)
            Idj = Id(J)

            Dx = Rx(I)-Rx(J)
            Dy = Ry(I)-Ry(J)
            Dz = Rz(I)-Rz(J)
            
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

            If(R2.Lt.Rcut2(Idi,Idj)) Then
               R2   = Sig2(Idi,Idj)/R2
               R2   = R2*R2*R2
               Uold = Uold + Eps(Idi,Idj)*R2*(R2-1.0d0) - Ecut(Idi,Idj)
               
               Virold = Virold  + 12.0d0*Eps(Idi,Idj)*R2*(R2-0.5d0)
            Endif
         Enddo
      Enddo
      
      Return
      End
