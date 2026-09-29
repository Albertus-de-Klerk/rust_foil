C=====================================================================
C     qfoil-rs reference driver: tabulates the QFoil BL closure
C     routines (xblsys.f) over a grid of inputs that crosses every
C     branch boundary.  Output: one line per call,
C       ROUTINE  inputs...  outputs...     (ES26.17E3)
C     Build/run: tools/fdrivers/run_closures.sh
C=====================================================================
      PROGRAM CLOSUR
      IMPLICIT REAL (A-H,M,O-Z)
      INCLUDE 'BLPAR.INC'
      DIMENSION HKS(24), RTS(12), MSQS(3), THS(3), AMS(5)
      DATA HKS / 1.02, 1.05, 1.3, 1.8, 2.2, 2.6, 3.0, 3.49, 3.5, 3.7,
     &           3.9, 4.0, 4.2, 4.35, 4.5, 5.0, 5.5, 5.8, 6.0, 7.0,
     &           8.0, 12.0, 15.0, 20.0 /
      DATA RTS / 20.0, 150.0, 199.0, 200.0, 250.0, 400.0, 401.0,
     &           1000.0, 3000.0, 1.0E4, 1.0E5, 1.0E6 /
      DATA MSQS / 0.0, 0.05, 0.3 /
      DATA THS / 1.0E-4, 1.0E-3, 1.0E-2 /
      DATA AMS / 0.0, 3.0, 8.2, 8.9, 9.5 /
C
      CALL BLPINI
C
      OPEN(10,FILE='closures.txt',STATUS='UNKNOWN')
  900 FORMAT(A8, 20ES26.17E3)
C
      DO 10 I=1, 24
        HK = HKS(I)
        DO 20 K=1, 3
          MSQ = MSQS(K)
          CALL HKIN(HK, MSQ, HKO, HKO_H, HKO_M)
          WRITE(10,900) 'HKIN', HK, MSQ, HKO, HKO_H, HKO_M
          CALL HCT(HK, MSQ, HC, HC_HK, HC_MSQ)
          WRITE(10,900) 'HCT', HK, MSQ, HC, HC_HK, HC_MSQ
   20   CONTINUE
        DO 30 J=1, 12
          RT = RTS(J)
          CALL DIL(HK, RT, DI, DI_HK, DI_RT)
          WRITE(10,900) 'DIL', HK, RT, DI, DI_HK, DI_RT
          CALL DILW(HK, RT, DI, DI_HK, DI_RT)
          WRITE(10,900) 'DILW', HK, RT, DI, DI_HK, DI_RT
          DO 40 K=1, 3
            MSQ = MSQS(K)
            CALL HSL(HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ)
            WRITE(10,900) 'HSL', HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ
            CALL HST(HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ)
            WRITE(10,900) 'HST', HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ
            CALL CFL(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
            WRITE(10,900) 'CFL', HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ
            CALL CFT(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
            WRITE(10,900) 'CFT', HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ
   40     CONTINUE
          DO 50 L=1, 3
            TH = THS(L)
            CALL DAMPL(HK, TH, RT, AX, AX_HK, AX_TH, AX_RT)
            WRITE(10,900) 'DAMPL', HK, TH, RT, AX, AX_HK, AX_TH, AX_RT
            CALL DAMPL2(HK, TH, RT, AX, AX_HK, AX_TH, AX_RT)
            WRITE(10,900) 'DAMPL2', HK, TH, RT, AX, AX_HK, AX_TH, AX_RT
   50     CONTINUE
   30   CONTINUE
   10 CONTINUE
C
C---- DIT over a small grid
      DO 60 I=1, 3
        HS = 1.4 + 0.3*FLOAT(I)
        DO 60 J=1, 3
          US = 0.3*FLOAT(J)
          DO 60 K=1, 3
            CF = 0.001*FLOAT(K)
            ST = 0.05*FLOAT(K+J)
            CALL DIT(HS, US, CF, ST, DI, DI_HS, DI_US, DI_CF, DI_ST)
            WRITE(10,900) 'DIT', HS, US, CF, ST,
     &                    DI, DI_HS, DI_US, DI_CF, DI_ST
   60 CONTINUE
C
C---- AXSET: pairs of stations, both amplification models
      DO 70 IDAMP=0, 1
      DO 70 I=2, 24, 3
        HK1 = HKS(I)
        HK2 = HKS(MIN(I+1,24))
        DO 70 J=4, 12, 2
          RT1 = RTS(J)
          RT2 = RTS(MIN(J+1,12))
          DO 70 L=1, 5
            A1 = AMS(L)
            A2 = AMS(MIN(L+1,5))
            T1 = 1.0E-3
            T2 = 1.3E-3
            CALL AXSET(HK1, T1, RT1, A1, HK2, T2, RT2, A2, 9.0, IDAMP,
     &                 AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &                     AX_HK2, AX_T2, AX_RT2, AX_A2)
            WRITE(10,900) 'AXSET', FLOAT(IDAMP), HK1, T1, RT1, A1,
     &                    HK2, T2, RT2, A2, AX, AX_HK1, AX_T1, AX_RT1,
     &                    AX_A1, AX_HK2, AX_T2, AX_RT2, AX_A2
   70 CONTINUE
C
      CLOSE(10)
      END
