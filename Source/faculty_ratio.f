      Function Faculty_Ratio(I,J)
      Implicit None

      Integer I,J,K
      Double Precision Faculty_Ratio
      
      Faculty_Ratio = 1.0d0

      If(I.Lt.0.Or.J.Lt.0) Stop "Error Faculty_Ratio"
      
      If(I.Gt.J) Then
         Do K=(J+1),I
            Faculty_Ratio = Faculty_Ratio*Dble(K)
         Enddo
      Elseif(J.Gt.I) Then
         Do K=(I+1),J
            Faculty_Ratio = Faculty_Ratio/Dble(K)
         Enddo
      Endif

      Return
      End
