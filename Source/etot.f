      Subroutine Etot(Ib,Vir,Upot)
      Implicit None

      Include 'commons.inc'

C     Compute The Energy Total Of Box Ib

      Integer I,J,Ib,Idi,Idj
      Double Precision Vir,Upot,Dx,Dy,Dz,R2,Bx,Hbx

      Upot = 0.0d0
      Vir  = 0.0d0

      Bx  = Box(Ib)
      Hbx = 0.5d0*Box(Ib)

      Do I=1,Npart-1

         If(Ibox(I).Eq.Ib) Then

            Idi = Id(I)

            Do J=I+1,Npart

               If(Ibox(J).Eq.Ib) Then

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
                     Upot = Upot + Eps(Idi,Idj)*R2*(R2-1.0d0) -
     &                    Ecut(Idi,Idj)
                     
                     Vir  = Vir + 12.0d0*Eps(Idi,Idj)*R2*(R2-0.5d0)
                  Endif
               Endif
            Enddo
         Endif
      Enddo

      Return
      End
