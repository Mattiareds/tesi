c
      SUBROUTINE RMARIN(IJ,KL)
      COMMON/RASET1/U(97),C,CD,CM,I97,J97
      I=MOD(IJ/177,177)+2
      J=MOD(IJ,177)+2
      K=MOD(KL/169,178)+1
      L=MOD(KL,169)
      DO II=1,97
      X=0.
      Y=.5
        DO  JJ=1,24
          M=MOD(MOD(I*J,179)*K,179)
          I=J
          J=K
          K=M
          L=MOD(53*L+1,169)
          IF (MOD(L*M,64).GE.32) X=X+Y
          Y=.5*Y
        enddo
        U(II)=X
      enddo
      C=362436./16777216.
      CD=7654321./16777216.
      CM=16777213./16777216.
      I97=97
      J97=33
      RETURN
      END
      SUBROUTINE RANMAR(zrvec)
      Real*8 zrvec
      COMMON/RASET1/U(97),C,CD,CM,I97,J97
      len=1
      DO IVEC=1,LEN
         UNI=U(I97)-U(J97)
         IF (UNI.LT.0.) UNI=UNI+1.
         U(I97)=UNI
         I97=I97-1
         IF (I97.EQ.0) I97=97
         J97=J97-1
         IF (J97.EQ.0) J97=97
         C=C-CD
         IF (C.LT.0.) C=C+CM
         UNI=UNI-C
         IF (UNI.LT.0.) UNI=UNI+1.
         zrvec=UNI
      enddo
      RETURN
      END
