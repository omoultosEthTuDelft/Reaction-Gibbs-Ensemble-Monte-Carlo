      Program Mc
      Implicit None

C     Monte Carlo Gibbs Ensemble With Chemical Reaction

      Include 'commons.inc'

      Logical Linit,Lnpt
      Integer Ncycle,Ninit,Nmove,I,J,Jj,Icycle,Icycle2
      Double Precision Pdisp,Pswap,Pvol,Preact,Temp,Deltax,Deltav,V,E,Rm
     $     ,Ran_Uniform,Dummy,Avv1,Avv2,Avs1,Avs2,Avd1,Avd2
     $     ,Avr1(Maxreact),Avr2(Maxreact),Epsilon(Maxcomp)
     $     ,Sigma(Maxcomp),R2

      Avv1 = 0.0d0
      Avv2 = 0.0d0
      Avs1 = 0.0d0
      Avs2 = 0.0d0
      Avd1 = 0.0d0
      Avd2 = 0.0d0

      Do I=1,Maxreact
         Avr1(I) = 0.0d0
         Avr2(I) = 0.0d0
      Enddo

      Rm = Ran_Uniform()
      
      Open(21,File="Input",Status="Unknown")
      Read(21,*)
      Read(21,*) Ncycle,Ninit,Linit,Temp
      Read(21,*)
      Read(21,*) Pdisp,Pswap,Pvol,Preact
      Read(21,*)
      Read(21,*) Deltax,Deltav,Lnpt,Press
      Read(21,*)
      Read(21,*) Box(1),Box(2)
      Read(21,*)
      Read(21,*) Ncomp

      If(Lnpt) Then
         Pswap = 0.0d0
         Nibox = 1
      Else
         Press = 0.0d0
         Nibox = 2
      Endif
      
      If(Press.Lt.0.0d0) Press = 0.0d0
      
      If(Ncomp.Lt.1.Or.Ncomp.Gt.Maxcomp) Stop "Error Maxcomp"
      
      Read(21,*)
      Read(21,*) (Nid(I,1),I=1,Ncomp)
      Read(21,*) (Nid(I,2),I=1,Ncomp)
      Read(21,*)
      Read(21,*) (Epsilon(I),I=1,Ncomp)
      Read(21,*) (Sigma(I),I=1,Ncomp)
      Read(21,*) (lnQ(I),I=1,Ncomp)
      Read(21,*)
      Read(21,*) Nreact

      If(Lnpt) Then
         Do I=1,Ncomp
            Nid(I,2) = 0
         Enddo
      Endif
      
      If(Ncomp.Eq.1) Nreact = 0
      
      If(Nreact.Lt.0.Or.Nreact.Gt.Maxreact)
     &     Stop "Error Nreact"

      Do J=1,Nreact
         Read(21,*) (Nstoi(Jj,J),Jj=1,Ncomp)
      Enddo
      Close(21)

      If(Temp.Le.0.0d0)   Stop "Error Temperature !!!"
      If(Pdisp.Le.0.0d0)  Pdisp  = 0.0d0
      If(Pvol.Le.0.0d0)   Pvol   = 0.0d0
      If(Pswap.Le.0.0d0)  Pswap  = 0.0d0
      If(Preact.Le.0.0d0) Preact = 0.0d0
      If(Deltax.Lt.0.0d0) Stop "Error Deltax !!"
      If(Deltav.Lt.0.0d0) Stop "Error Deltav !!"
      If(Ncycle.Lt.100)   Stop "Minimal 100 Cycles !!!"

      Beta   = 1.0d0/Temp
      Rcutsq = 0.0d0
      
      Do J=1,Ncomp
         Do Jj=1,Ncomp
            Sig2(J,Jj)  = (0.5d0*(Sigma(J)+Sigma(Jj)))**2
            Eps(J,Jj)   = 4.0d0*Dsqrt(Abs(Epsilon(J)*Epsilon(Jj)))
            Rcut2(J,Jj) = (2.5d0*Dsqrt(Sig2(J,Jj)))**2

            R2 = (Sig2(J,Jj)/Rcut2(J,jj))**3
            
            Ecut(J,Jj) = Eps(J,Jj)*R2*(R2-1.0d0)

            Rcutsq = Max(Rcutsq,Rcut2(J,Jj))
         Enddo
      Enddo
      Write(6,*) 'Units printed here are correct if epsilon/kB is in [K], sigma is in [Angstrom] and T is in [K].'
      Write(6,*)
      Write(6,*) 'Ncycle               : ',Ncycle
      Write(6,*) 'Ninit                : ',Ninit
      Write(6,*) 'Linit                : ',Linit
      Write(6,*) 'Temp                 : ',Temp, ' [K]'
      Write(6,*) 'Deltax               : ',Deltax, ' [Angstrom]'
      Write(6,*) 'Deltav               : ',Deltav, ' [cubic Angstrom]'
      Write(6,*) 'Lnpt                 : ',Lnpt
      Write(6,*) 'Nibox                : ',Nibox
      Write(6,*) 'Pressure             : ',Press, ' [kPa]'

      Press = Press / to_kPa

      Write(6,*)
      Write(6,*) 'Ncomp                : ',Ncomp
      Write(6,*)

      Do J=1,Ncomp
         Do Jj=1,Ncomp
            Write(6,*) 'Epsilon/kB    : ',J,Jj,Eps(J,Jj)*0.25d0, ' [K]'
            Write(6,*) 'Sigma         : ',J,Jj,Dsqrt(Sig2(J,Jj)), ' [Angstrom]'
            Write(6,*) 'Rcut          : ',J,Jj,Dsqrt(Rcut2(J,Jj)), ' [Angstrom]'
            Write(6,*) 'Ecut          : ',J,Jj,Rgas * Ecut(J,Jj), ' [kJ/mol]'
         Enddo
      Enddo

      Write(6,*)
      Write(6,*) 'Rcutmax              : ',Dsqrt(Rcutsq), ' [Angstrom]'
      Write(6,*)
      Write(6,*)
      Write(6,*) 'ln(Q)                : ',(lnQ(J),J=1,Ncomp)
      Write(6,*)
      Write(6,*) 'Nreact               : ',Nreact
      Write(6,*)

      Do J=1,Nreact
         Write(6,*) 'Reaction    : ',J
         Write(6,*) 'Nstoi       : ',(Nstoi(Jj,J),Jj=1,Ncomp)
      Enddo

      Write(6,*)
      Write(6,*)
            
      If(Linit) Then
         Write(6,*)
         Write(6,*) 'Generate Initial Coordinates'
         Write(6,*)

         If(Dsqrt(Rcutsq).Gt.0.5d0*Min(Box(1),Box(2)))
     &        Stop "Boxes Too Small !!!"
         
         Call Init
      Else
         Write(6,*)
         Write(6,*) 'Read Coordinates From Disk'
         Write(6,*)

         Open(21,File="Coordold",Status="Unknown")
         Read(21,*) Box(1),Box(2)
         Read(21,*) Npart,Ncomp

         If(Npart.Le.0.Or.Npart.Gt.Maxpart.Or.
     &        Ncomp.Le.0.Or.Ncomp.Gt.Maxcomp) Then

            Write(6,*) Maxpart,Npart
            Write(6,*) Maxcomp,Ncomp
            Stop
         Endif
         
         Read(21,*) Npbox(1),Npbox(2)
         Read(21,*) (Nid(J,1),J=1,Ncomp)
         Read(21,*) (Nid(J,2),J=1,Ncomp)

         If(Dsqrt(Rcutsq).Gt.0.5d0*Min(Box(1),Box(2)))
     &        Stop "Boxes Too Small !!!"

         Do I=1,Npart
            Read(21,*) Rx(I),Ry(I),Rz(I),Ibox(I),Id(I)
         Enddo
      Endif

      Call Check
      
      Write(6,*) 'Box 1                : ',Box(1), ' [Angstrom]'
      Write(6,*) 'Box 2                : ',Box(2), ' [Angstrom]'
      Write(6,*) 'Npart                : ',Npart
      Write(6,*) 'Npbox 1              : ',Npbox(1)
      Write(6,*) 'Npbox 2              : ',Npbox(2)
      Write(6,*) '<Rho1>               : ',Dble(Npbox(1))/(Box(1)**3), ' [particles/cubic Angstrom]'
      Write(6,*) '<Rho2>               : ',Dble(Npbox(2))/(Box(2)**3), ' [particles/cubic Angstrom]'
      Write(6,*)
      Write(6,*)
      
      Do Jj=1,Ncomp
         Write(6,*) 'Nid Box1    : ',Jj,Nid(Jj,1)
         Write(6,*) 'Rho Box1    : ',Jj,Dble(Nid(Jj,1))/(Box(1))**3, ' [particles/cubic Angstrom]'
      Enddo

      Write(6,*)
      Write(6,*)

      Do Jj=1,Ncomp
         Write(6,*) 'Nid Box2    : ',Jj,Nid(Jj,2)
         Write(6,*) 'Rho Box2    : ',Jj,Dble(Nid(Jj,2))/(Box(2))**3, ' [particles/cubic Angstrom]'
      Enddo

      Write(6,*)
      Write(6,*)

      Dummy  = Pdisp + Pswap + Pvol + Preact
      Pdisp  = Pdisp/Dummy
      Pswap  = Pswap/Dummy
      Pvol   = Pvol/Dummy
      Preact = Preact/Dummy

      Write(6,*) 'Pdisp                : ',Pdisp
      Write(6,*) 'Pswap                : ',Pswap
      Write(6,*) 'Pvol                 : ',Pvol
      Write(6,*) 'Preact               : ',Preact
      Write(6,*)
     
      Etotal(1) = 0.0d0
      Etotal(2) = 0.0d0
      Vtotal(1) = 0.0d0
      Vtotal(2) = 0.0d0

      Do I=1,2
         Call Etot(I,V,E)

         Vtotal(I) = V
         Etotal(I) = E
      Enddo

      Write(6,*)
      Write(6,*) 'Initial Energy Box 1 : ', Rgas * Etotal(1), ' [kJ/mol]'
      Write(6,*) 'Initial Energy Box 2 : ', Rgas * Etotal(2), ' [kJ/mol]'
      Write(6,*)
      Write(6,*) 'Initial Virial Box 1 : ', to_kPa * Vtotal(1), ' [kPa]'
      Write(6,*) 'Initial Virial Box 2 : ', to_kPa * Vtotal(2), ' [kPa]'
      Write(6,*)

C     Start Of The Simulation

      Write(6,*)
      Write(6,*)
      Write(6,*) 'The Simulation Is Running.....'
      Write(6,*)
      Write(6,*)

      Call Sample(1)

      Open(22,File="Traject.xyz",Status="Unknown")
      Open(23,File="Results",Status="Unknown")
    
      Do Icycle=1,Ncycle

         Nmove = Max(20,Npart)

         Do Icycle2=1,Nmove

C     Select Trial Move At Random

            Rm = Ran_Uniform()

            If(Rm.Lt.Pvol) Then

               If(Lnpt) Then
                  Call Volume_NPT(Avv1,Avv2,Deltav)
               Else
                  Call Volume_Gibbs(Avv1,Avv2,Deltav)
               Endif

            Elseif(Rm.Lt.Pswap+Pvol) Then

               Call Swap_Gibbs(Avs1,Avs2)
                          
            Elseif(Rm.Lt.Pswap+Pvol+Pdisp) Then

               Call Move(Avd1,Avd2,Deltax)

            Else

               Call Reaction(Avr1,Avr2)
               
            Endif

            If(Icycle.Gt.Ninit) Call Sample(2)

         Enddo

         If(Mod(Icycle,100).Eq.0) Then
            Write(23,'(7e20.10)')
     &           Dble(Icycle),
     &           Dble(Npbox(1))/(Box(1)**3),
     &           Dble(Npbox(2))/(Box(2)**3),
     &           Etotal(1)/Max(0.5d0,Dble(Npbox(1))),
     &           Etotal(2)/Max(0.5d0,Dble(Npbox(2))),
     &           to_kPa*(Dble(Npbox(1))/(Beta*(Box(1)**3)) + Vtotal(1)/(3.0d0*(Box(1)**3))),
     &           to_kPa*(Dble(Npbox(2))/(Beta*(Box(2)**3)) + Vtotal(2)/(3.0d0*(Box(2)**3)))
         Endif
      
         If(Icycle.Gt.Ninit) Then
            Call Sample(3)
         Endif

         If(Mod(Icycle,Ncycle/100).Eq.0) Then
            Write(22,*) Npart
            Write(22,*)

            Do J=1,Npart
               If(Ibox(J).Eq.2) Then
                  Dummy = 4.0d0*(Box(1) + 2.0d0)
               Else
                  Dummy = 0.0d0
               Endif
               
               Write(22,'(A,3f15.5)') 'Ar  ',4.0d0*Rx(J)+Dummy,
     &              4.0d0*Ry(J),4.0d0*Rz(J)
            Enddo 
         Endif 
      Enddo

      Close(22)
      Close(23)
   
      Call Sample(4)
     
      Write(6,*) 'Frac. Acc. Displ.    : ',Avd1/Max(0.5d0,Avd2)
      Write(6,*) 'Frac. Acc. Swap      : ',Avs1/Max(0.5d0,Avs2)
      Write(6,*) 'Frac. Acc. Volume    : ',Avv1/Max(0.5d0,Avv2)
           
      Do I=1,Nreact
         Write(6,*) 'Frac. Acc. React     : ',Avr1(I)/Max(0.5d0,Avr2(I))
      Enddo
      
      Write(6,*)

      Open(21,File="Coordnew",Status="Unknown")
      Write(21,*) Box(1),Box(2)
      Write(21,*) Npart,Ncomp
      Write(21,*) Npbox(1),Npbox(2)
      Write(21,*) (Nid(J,1),J=1,Ncomp)
      Write(21,*) (Nid(J,2),J=1,Ncomp)
      
      Do I=1,Npart
         Write(21,'(3e20.10,2i5)') Rx(I),Ry(I),Rz(I),Ibox(I),Id(I)
      Enddo
      Close(21)

C     Check Energy Calculation

      Do I=1,2
         Call Etot(I,V,E)
       
         Write(6,*)
         Write(6,*) 'Box                  : ',I
         Write(6,*) 'Particles            : ',Npbox(I)
         Write(6,*) 'Box                  : ',Box(I), ' [Angstrom]'
         Write(6,*) 'Volume               : ',Box(I)**3, ' [cubic Angstrom]'
         Write(6,*) 'Rho                  : ',Dble(Npbox(I))/(Box(I)**3), ' [particles/cubic Angstrom]'
         Write(6,*) 'Energy               : ',Rgas * E, ' [kJ/mol]'
         Write(6,*) 'Energy Sim.          : ',Rgas * Etotal(I), ' [kJ/mol]'
         Write(6,*) 'Diff #####           : ',Rgas * Dabs(Etotal(I)-E), ' [kJ/mol]'
         Write(6,*) 'Virial               : ',to_kPa * V, ' [kPa]'
         Write(6,*) 'Virial Sim.          : ',to_kPa * Vtotal(I), ' [kPa]'
         Write(6,*) 'Diff #####           : ',to_kPa * Dabs(Vtotal(I)-V), ' [kPa]'
         Write(6,*)
      Enddo

      Write(6,*)
      Write(6,*)
      
      Do Jj=1,Ncomp
         Write(6,*) 'Nid Box1 : ',Jj,Nid(Jj,1)
         Write(6,*) 'Rho Box1 : ',Jj,Dble(Nid(Jj,1))/(Box(1))**3
     &        ,' [particles/cubic Angstrom]'
      Enddo

      Write(6,*)
      Write(6,*)

      Do Jj=1,Ncomp
         Write(6,*) 'Nid Box2 : ',Jj,Nid(Jj,2)
         Write(6,*) 'Rho Box2 : ',Jj,Dble(Nid(Jj,2))/(Box(2))**3
     &        ,' [particles/cubic Angstrom]'
      Enddo

      Call Check
      
      Write(6,*)
      Write(6,*)

      Stop
      End
