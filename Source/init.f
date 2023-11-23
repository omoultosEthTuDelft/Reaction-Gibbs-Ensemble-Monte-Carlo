      Subroutine Init
      Implicit None

      Include 'commons.inc'

C     Generate An Initial Configuration

      Logical Laccept
      Integer I,J,Jj,K,Ipart,Ib,Idp,Selectint
      Double Precision Ran_Uniform,Xi,Yi,Zi,Rxn,Ryn,Rzn,Unew,Uold
     $     ,Virnew,Virold

      I = 0

      Npbox(1) = 0
      Npbox(2) = 0
      
      Do J=1,2
         Do Jj=1,Ncomp

            If(Nid(Jj,j).Lt.0) Stop "Error Nid Init"
            
            Do K=1,Nid(Jj,J)
     
               I = I + 1

               Npart = I

               If(Npart.Gt.Maxpart) Stop "Increase Maxpart"
               
               Id(I)   = Jj
               Ibox(I) = J

               Npbox(J) = Npbox(J) + 1
                           
               Rx(I) = Box(J)*Ran_Uniform()
               Ry(I) = Box(J)*Ran_Uniform()
               Rz(I) = Box(J)*Ran_Uniform()
            Enddo
         Enddo
      Enddo

      If(Npart.Eq.0) Stop "Error Maxpart"
      
C     Monte Carlo Displacements To Remove Initial Overlaps

      Do J=1,Npart*50
         Ipart = Selectint(Npart)
         Ib    = Ibox(Ipart)
         Idp   = Id(Ipart)

         Rxn = Rx(Ipart) + Ran_Uniform() - 0.5d0
         Ryn = Ry(Ipart) + Ran_Uniform() - 0.5d0
         Rzn = Rz(Ipart) + Ran_Uniform() - 0.5d0

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

         Call Epart(Ib,Virold,Uold,Xi,Yi,Zi,Ipart,Idp)
         Call Epart(Ib,Virnew,Unew,Rxn,Ryn,Rzn,Ipart,Idp)

         Call Accept(Dexp(-0.5d0*(Unew-Uold)),Laccept)

         If(Laccept) Then
            Rx(Ipart) = Rxn
            Ry(Ipart) = Ryn
            Rz(Ipart) = Rzn
         Endif
      Enddo

      Return
      End
