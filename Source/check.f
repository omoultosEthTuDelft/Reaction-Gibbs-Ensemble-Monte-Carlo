      Subroutine Check
      Implicit None

      Include 'commons.inc'

      Integer Dummy,I,J,Jj
      
      If(Npart.Ne.(Npbox(1)
     $     +Npbox(2)).Or.Npart.Le.0.Or.Npart.Gt.Maxpart
     &     .Or.Npbox(1).Lt.0.Or.Npbox(2).Lt.0) Then
         
         Write(6,*) 'Error 1'
         Write(6,*) Maxpart,Npart,Npbox(1),Npbox(2)
         Stop
      Endif

      Do Jj=1,2
         Dummy = 0
         Do J=1,Ncomp
            If(Nid(J,Jj).Lt.0) Then
               Write(6,*) 'Error 2'
               Write(6,*) Jj,J,Nid(J,Jj)
               Stop
            Endif
            
            Dummy = Dummy + Nid(J,Jj)
         Enddo

         If(Dummy.Ne.Npbox(Jj)) Then
            Write(6,*) 'Error 3'
            Write(6,*) Jj,Npbox(Jj),Dummy
            Stop
         Endif
      Enddo

      Do Jj=1,2

         Do J=1,Ncomp

            Dummy = 0
            
            Do I=1,Npart
               If(Id(I).Le.0.Or.Id(I).Gt.Ncomp) Then
                  Write(6,*) 'Error 4'
                  Write(6,*) Jj,J,I,Id(I)
                  Stop
               Endif
               
               If(Id(I).Eq.J.And.Ibox(I).Eq.Jj)
     &              Dummy = Dummy + 1
            Enddo

            If(Dummy.Ne.Nid(J,Jj)) Then
               Write(6,*) 'Error 5'
               Write(6,*) Jj,J,Dummy,Nid(J,Jj)
               Stop
            Endif
         Enddo
      Enddo

      Do Jj=1,2
         Dummy = 0
            
         Do I=1,Npart

            If(Ibox(I).Ne.1.And.Ibox(I).Ne.2) Then
               Write(6,*) 'Error 6'
               Write(6,*) I,Ibox(I)
               Stop
            Endif

            If(Ibox(I).Eq.Jj) Dummy = Dummy + 1
         Enddo

         If(Dummy.Ne.Npbox(Jj)) Then
            Write(6,*) 'Error 7'
            Write(6,*) Jj,Dummy,Npbox(Jj)
            Stop
         Endif
      Enddo
            
      Return
      End
      
