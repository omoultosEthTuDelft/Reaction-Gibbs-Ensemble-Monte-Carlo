      Subroutine Energy_Add(Ib,Unew,Virnew,Nadd,Idadd,
     &     Rxadd,Ryadd,Rzadd,Nremove,Iremove)
      
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Of Particles That Need To Be Removed
      
      Integer Ib,Nremove,Iremove(Maxpart),I,J,Idi,Idj,Nadd
     $     ,Idadd(Maxpart)

      Double Precision Unew,Virnew,R2,Dx,Dy,Dz,Bx,Hbx,Rxadd(Maxpart)
     $     ,Ryadd(Maxpart),Rzadd(Maxpart)

      Logical Lok(Maxpart)

C     Check For The Correct Box
      
      Do I=1,Npart
         If(Ibox(I).Eq.Ib) Then
            Lok(I) = .True.
         Else
            Lok(I) = .False.
         Endif
      Enddo

C     Exclude Interactions With The Reactants Than Will Be Removed

      Do I=1,Nremove
         Lok(Iremove(I)) = .False.
      Enddo

      Unew   = 0.0d0
      Virnew = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

C     For Each Particle That Need To Be Added, Calculate The Interactions With
C     All Other Particles, Except The Ones That Need To Be Removed
      
      Do I=1,Npart

         If(Lok(I)) Then

            Idi = Id(I)

            Do J=1,Nadd

               Idj = Idadd(J)
                              
               Dx = Rx(I)-Rxadd(J)
               Dy = Ry(I)-Ryadd(J)
               Dz = Rz(I)-Rzadd(J)
            
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
                  Unew = Unew + Eps(Idi,Idj)*R2*(R2-1.0d0) - Ecut(Idi,Idj)
                  
                  Virnew = Virnew  + 12.0d0*Eps(Idi,Idj)*R2*(R2-0.5d0)
               Endif
            Enddo
         Endif
      Enddo

C     Compute The Energy Between Particles That Need To Be Added
      
      Do I=1,Nadd-1
        
         Idi = Idadd(I)

         Do J=I+1,Nadd

            Idj = Idadd(J)

            Dx = Rxadd(I)-Rxadd(J)
            Dy = Ryadd(I)-Ryadd(J)
            Dz = Rzadd(I)-Rzadd(J)
            
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
               Unew = Unew + Eps(Idi,Idj)*R2*(R2-1.0d0) - Ecut(Idi,Idj)
               
               Virnew = Virnew  + 12.0d0*Eps(Idi,Idj)*R2*(R2-0.5d0)
            Endif
         Enddo
      Enddo
      
      Return
      End
