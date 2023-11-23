      Subroutine Sample(Ic)
      Implicit None

C     Sample Averages

      Include 'commons.inc'

      Integer I,J,Jj,Ic
      Double Precision Av1(5,2),Av2(Maxcomp,2),Av3(Maxcomp,2),Av4
     $     ,Avmu1(Maxcomp,2),Avmu2(Maxcomp,2),Ran_Uniform,Xi,Yi,Zi,Up
     $     ,Dummy,Eq

      Save Av1,Av2,Av3,Av4,Avmu1,Avmu2

      If(Ic.Eq.1) Then

C     Set Counters To Zero

         Av4 = 0.0d0
     
         Do I=1,5
            Do J=1,2
               Av1(I,J) = 0.0d0
            Enddo
         Enddo

         Do Jj=1,Ncomp
            Do J=1,2
               Av2(Jj,J) = 0.0d0
               Av3(Jj,J) = 0.0d0

               Avmu1(Jj,J) = 0.0d0
               Avmu2(Jj,J) = 0.0d0
            Enddo
         Enddo
        
      Elseif(Ic.Eq.2) Then

C     Sample Ensemble Averages

         Av4 = Av4 + 1.0d0

         Do J=1,Nibox
            Av1(1,J) = Av1(1,J) + Etotal(J)
            
            Av1(2,J) = Av1(2,J) + Dble(Npbox(J))/(Beta*(Box(J)**3)) + 
     &           Vtotal(J)/(3.0d0*(Box(J)**3))
            
            Av1(3,J) = Av1(3,J) + Dble(Npbox(J))/(Box(J)**3)
            
            Av1(4,J) = Av1(4,J) + Box(J)**3

            Av1(5,J) = Av1(5,J) + Dble(Npbox(J))

            Do Jj=1,Ncomp
               Av2(Jj,J) = Av2(Jj,J) + Dble(Nid(Jj,J))
               Av3(Jj,J) = Av3(Jj,J) + Dble(Nid(Jj,J))/(Box(J)**3)
            Enddo
         Enddo

      Elseif(Ic.Eq.3) Then

C     Sample Chemical Potential

         Do J=1,Nibox
            Do Jj=1,Ncomp
               Do I=1,Npart
                  Xi = Ran_Uniform()*Box(J)
                  Yi = Ran_Uniform()*Box(J)
                  Zi = Ran_Uniform()*Box(J)
               
                  Call Epart(J,Dummy,Up,Xi,Yi,Zi,0,Jj)

                  Avmu1(Jj,J) = Avmu1(Jj,J) + Dexp(-Beta*Up)
                  Avmu2(Jj,J) = Avmu2(Jj,J) + 1.0d0
               Enddo
            Enddo
         Enddo
         
      Else

C     Write Averages

         Do Jj=1,Nibox
            Write(6,*)
            Write(6,*)
            Write(6,*) '##############################################################################'
            Write(6,*)
            Write(6,*)
            Write(6,*) 'Averages Box ',Jj
            Write(6,*)
            Write(6,*) '<E>                  : ',Rgas * Av1(1,Jj)/Av4, ' [kJ/mol]'
            Write(6,*) '<P>                  : ',to_kPa * Av1(2,Jj)/Av4, ' [kPa]'
            Write(6,*) '<Rho>                : ',Av1(3,Jj)/Av4, ' [particles/cubic Angstrom]'
            Write(6,*) '<V>                  : ',Av1(4,Jj)/Av4, ' [cubic Angstrom]'
            Write(6,*) '<N>                  : ',Av1(5,Jj)/Av4
            Write(6,*)

            Do J=1,Ncomp
               Write(6,*) 'Component            : ',J
               Write(6,*) 'N                    : ',Av2(J,Jj)/Av4
               Write(6,*) 'Rho                  : ',Av3(J,Jj)/Av4, ' [particles/cubic Angstrom]'
            
               Write(6,*) 'Mu                   : ',
     &              Rgas * (Dlog(Av3(J,Jj)/Av4)/Beta -
     &              Dlog(Avmu1(J,Jj)/Avmu2(J,Jj))/Beta -
     &              lnQ(J)/Beta), ' [kJ/mol]'
            
               Write(6,*) 'Mu_Ex                : ',
     &              Rgas * (Dlog(Av3(J,Jj)/Av4)/Beta),
     &              ' [kJ/mol]'

               Write(6,*) 'Mu_Ig                : ',
     &              Rgas * 
     &              (-Dlog(Avmu1(J,Jj)/Avmu2(J,Jj))/Beta),
     &              ' [kJ/mol]'

               Write(6,*) 'Mu_Intra             : ',
     &              Rgas * (-lnQ(J)/Beta), ' [kJ/mol]'
            Enddo

            Write(6,*)
            Write(6,*) '##############################################################################'
            Write(6,*)
            Write(6,*)
         Enddo

C     Write Actual and Ideal Reaction Equilibrium Constant

         Do I=1,Nreact
            Do Jj=1,Nibox
              
               Write(6,*)
               Write(6,*)
               Write(6,*) 'Reaction             : ',I
               Write(6,*) 'Box                  : ',Jj

               Eq = 0.0d0

               Do J=1,Ncomp
                  Eq = Eq + Dble(Nstoi(J,I))*lnQ(J)
               Enddo

               Eq = Dexp(Eq)

               Write(6,*) 'Ideal K              : ',Eq
              
               Eq = 1.0d0

               Do J=1,Ncomp
                  Eq = Eq*Dexp(Dble(Nstoi(J,I))*Dlog(Av3(J,Jj)/Av4))
               Enddo

               Write(6,*) 'Real K               : ',Eq
               Write(6,*)
            Enddo
         Enddo
      Endif

      Return
      End
