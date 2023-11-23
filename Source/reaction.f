      Subroutine Reaction(Avr1,Avr2)
      Implicit None

      Include 'commons.inc'
      
      Double Precision Avr1(Maxreact),Avr2(Maxreact),Ran_Uniform
     $     ,Faculty_Ratio,Uold,Virold,Unew,Virnew,Rxadd(Maxpart)
     $     ,Ryadd(Maxpart),Rzadd(Maxpart),Factor,Rxold(Maxpart)
     $     ,Ryold(Maxpart),Rzold(Maxpart),SumQ

      Integer I,Ii,J,Is,Ib,Selectint,Ireact,Nadd,Nremove
     $     ,Iremove(Maxpart),Idadd(Maxpart),Ipart,Idold(Maxpart)
     $     ,Iboxold(Maxpart),Npartold

      Logical Laccept,Lremove(Maxpart)
      
      If(Nreact.Eq.0) Return

C     Select Direction Of The Reaction At Random
      
      If(Ran_Uniform().Lt.0.5d0) Then
         Is = 1
      Else
         Is = -1
      Endif

C     Select Box And Reaction At Random
      
      Ib     = Selectint(2)
      Ireact = Selectint(Nreact)

      Avr2(Ireact) = Avr2(Ireact) + 1.0d0
      
C     Select Reactant Molecules To Be Deleted

      Nremove = 0

      Do I=1,Maxpart
         Iremove(I) = 0
      Enddo
           
      Do I=1,Ncomp
         
         If(Is*Nstoi(I,Ireact).Lt.0) Then

C     Check If Molecules Are Available Otherwise Reject
            
            If(Nid(I,Ib).Lt.Abs(Nstoi(I,Ireact))) Return
                  
C     Select Random Molecules To Be Removed
            
            Do Ii=1,Abs(Nstoi(I,Ireact))

 1             Continue
               Ipart = Selectint(Npart)
               If(Id(Ipart).Ne.I) Goto 1
               If(Ibox(Ipart).Ne.Ib) Goto 1

C     Cannot Be A Particle That Was Already Selected
               
               Do J=1,Nremove
                  If(Iremove(J).Eq.Ipart) Goto 1
               Enddo

               Nremove          = Nremove + 1
               Iremove(Nremove) = Ipart
               
            Enddo
         Endif
      Enddo

C     Compute Energy Of Particles That Need To Be Removed

      Call Energy_Remove(Ib,Uold,Virold,Nremove,Iremove)

C     Add Particles At Random Positions In The Simulation Box

      Nadd = 0

      Do I=1,Ncomp

         If(Is*Nstoi(I,Ireact).Gt.0) Then
         
            Do Ii=1,Abs(Nstoi(I,Ireact))

               Nadd = Nadd + 1

               If(Nadd.Gt.Maxpart) Stop "Error Maxpart Nadd"

               Rxadd(Nadd) = Box(Ib)*Ran_Uniform()
               Ryadd(Nadd) = Box(Ib)*Ran_Uniform()
               Rzadd(Nadd) = Box(Ib)*Ran_Uniform()
               Idadd(Nadd) = I

            Enddo
         Endif
      Enddo

      Call Energy_Add(Ib,Unew,Virnew,Nadd,Idadd,Rxadd,Ryadd,Rzadd,
     &     Nremove,Iremove)

      Factor = Dexp(-Beta*(Unew-Uold))

      SumQ = 0.0d0
      
      Do I=1,Ncomp
         If(Nstoi(I,Ireact).Ne.0) Then
            Factor = Factor*Faculty_Ratio(Nid(I,Ib),
     &           (Nid(I,Ib)+(Is*Nstoi(I,Ireact))))
            
C                Dexp(Dble(Is)*Dble(Nstoi(I,Ireact))*
C                Dlog((Box(Ib)**3)*Q(I)))

            SumQ = SumQ + Dble(Is)*Dble(Nstoi(I,Ireact))*
     &           (3.0d0*Dlog(Box(Ib)) + lnQ(I))
         Endif
      Enddo

      Factor = Factor * Dexp(SumQ)
      
      Call Accept(Factor,Laccept)

      If(Laccept) Then
         Avr1(Ireact) = Avr1(Ireact) + 1.0d0

         Etotal(Ib) = Etotal(Ib) + Unew   - Uold
         Vtotal(Ib) = Vtotal(Ib) + Virnew - Virold

         Do I=1,Npart
            Lremove(I) = .False.
            Rxold(I)   = Rx(I)
            Ryold(I)   = Ry(I)
            Rzold(I)   = Rz(I)
            Idold(I)   = Id(I)
            Iboxold(I) = Ibox(I)
         Enddo

         Do I=1,Nremove
            Lremove(Iremove(I)) = .True.
         Enddo

         Npartold = Npart
         Npart    = 0
         Npbox(1) = 0
         Npbox(2) = 0

         Do I=1,Ncomp
            Nid(I,1) = 0
            Nid(I,2) = 0
         Enddo

         Do I=1,Npartold
            If(.Not.Lremove(I)) Then
               Npart = Npart + 1

               Rx(Npart)   = Rxold(I)
               Ry(Npart)   = Ryold(I)
               Rz(Npart)   = Rzold(I)
               Id(Npart)   = Idold(I)
               Ibox(Npart) = Iboxold(I)

               Npbox(Ibox(Npart))         = Npbox(Ibox(Npart))         + 1
               Nid(Id(Npart),Ibox(Npart)) = Nid(Id(Npart),Ibox(Npart)) + 1
            Endif
         Enddo
         
                           
         Do I=1,Nadd
            Npart = Npart + 1

            If(Npart.Gt.Maxpart) Stop "Error Maxpart React!"

            Npbox(Ib)        = Npbox(Ib)        + 1
            Nid(Idadd(I),Ib) = Nid(Idadd(I),Ib) + 1
            
            Rx(Npart) = Rxadd(I)
            Ry(Npart) = Ryadd(I)
            Rz(Npart) = Rzadd(I)

            Id(Npart)   = Idadd(I)
            Ibox(Npart) = Ib
         Enddo
      Endif
      
      Return
      End
      
