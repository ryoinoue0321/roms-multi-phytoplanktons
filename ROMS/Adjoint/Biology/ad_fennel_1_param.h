#include "fennel_bs1.h"

#if defined OXYGEN && defined DENITRIFICATION
!
!-----------------------------------------------------------------------
! Denitrification in anoxic water                       Okada 2014/02/13
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
            DO k=1,N(ng)
!             ==========================================================
              fac1=dtdays*DenitR(ng)
              cff1=MAX(Bio(i,k,iOxyg),0.0_r8)/K_Denit(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
# ifdef TDEPENDANCE
              fac2=DenitR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              cff3=cff2*fac1*fac2
# else
              cff3=cff2*fac1
# endif
!             ==========================================================
!!            Bio3(i,k,iNO3_)=Bio(i,k,iNO3_)
!>            tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio(i,k,iNO3_)*      &
!>   &                           tl_cff3)/(1.0_r8+cff3)
              adfac=ad_Bio(i,k,iNO3_)/(1.0_r8+cff3)
              ad_cff3=ad_cff3-Bio(i,k,iNO3_)*adfac
              ad_Bio(i,k,iNO3_)=adfac
# ifdef TDEPENDANCE
!>            tl_cff3=tl_cff2*fac1*fac2+                                &
!>   &                cff2*tl_fac1*fac2+                                &
!>   &                cff2*fac1*tl_fac2
              ad_fac2=ad_fac2+cff2*fac1*ad_cff3
              ad_fac1=ad_fac1+cff2*ad_cff3*fac2
              ad_cff2=ad_cff2+ad_cff3*fac1*fac2
              ad_cff3=0.0_r8

!>            tl_fac2=fac2*tl_Bio(i,k,itemp)*LOG(DenitR_t(ng))
#  ifndef UV_FIXED_TL
              ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+                      &
     &                          fac2*ad_fac2*LOG(DenitR_t(ng))
#  endif
              ad_fac2=0.0_r8
# else
!>            tl_cff3=tl_cff2*fac1+cff2*tl_fac1
              ad_fac1=ad_fac1+cff2*ad_cff3
              ad_cff2=ad_cff2+ad_cff3*fac1
              ad_cff3=0.0_r8
# endif
!>            tl_cff2=-cff2*cff2*tl_cff1
              ad_cff1=ad_cff1-cff2*cff2*ad_cff2
              ad_cff2=0.0_r8

!>            tlfac=(0.5_r8+SIGN(0.5_r8,Bio(i,k,iOxyg)))
!>            tl_cff1=(tlfac*tl_Bio(i,k,iOxyg)-cff1*tl_K_Denit)/        &
!>   &                K_Denit(ng)
              adfac=(0.5_r8+SIGN(0.5_r8,Bio(i,k,iOxyg)))
              ad_Bio(i,k,iOxyg)=ad_Bio(i,k,iOxyg)+                      &
     &                          adfac*ad_cff1/K_Denit(ng)
              ad_K_Denit=ad_K_Denit-cff1*ad_cff1/K_Denit(ng)
              ad_cff1=0.0_r8

!>            tl_fac1=dtdays*tl_DenitR
              ad_DenitR=ad_DenitR+dtdays*ad_fac1
              ad_fac1=0.0_r8
            END DO
          END DO
#endif
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            AttFac=0.0_r8
            IF (PARsur(i).gt.0.0_r8) THEN
!>            DO k=N(ng),1,-1
              DO k=1,N(ng)
!
! Compute the basic state PAR appropriate for each level.
!
!               ========================================================
                PAR=PARsur(i)
                DO kk=N(ng),k,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
#ifdef PHYT2
!                fac1=AttSW(ng)+AttChl(ng)*(Bio1(i,k,iChlo1)+Bio1(i,k,iChlo2))+AttFac
                  fac1=AttSW(ng)+AttChl(ng)*Bio1(i,kk,iChlo1)+AttFac
                  fac12=AttSW(ng)+AttChl(ng)*Bio1(i,kk,iChlo2)+AttFac

                  Att=fac1*(z_w(i,j,kk)-z_w(i,j,kk-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  fac2=Itop*(1.0_r8-ExpAtt)
                  PAR1=fac2/Att

                  Att2=fac12*(z_w(i,j,kk)-z_w(i,j,kk-1))
                  ExpAtt2=EXP(-Att2)
                  Itop2=PAR2
                  fac22=Itop2*(1.0_r8-ExpAtt2)
                  PAR12=fac22/Att2
                  PAR2=Itop2*ExpAtt2
#else
                  fac1=AttSW(ng)+AttChl(ng)*Bio1(i,kk,iChlo)+AttFac
                  Att=fac1*(z_w(i,j,kk)-z_w(i,j,kk-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  fac2=Itop*(1.0_r8-ExpAtt)
                  PAR1=fac2/Att
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
!               ========================================================
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
!>              tl_PAR=tl_Itop*ExpAtt+Itop*tl_ExpAtt
                ad_ExpAtt=ad_ExpAtt+Itop*ad_PAR
                ad_Itop=ad_Itop+ad_PAR*ExpAtt
                ad_PAR=0.0_r8

#ifdef PHYT2
!>              tl_PAR2=tl_Itop2*ExpAtt2+Itop2*tl_ExpAtt2
                ad_ExpAtt2=ad_ExpAtt2+Itop2*ad_PAR2
                ad_Itop2=ad_Itop2+ad_PAR2*ExpAtt2
                ad_PAR2=0.0_r8

#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
!               ========================================================
                fac=dtdays*NitriR(ng)
#ifdef OXYGEN
                fac2=MAX(Bio2(i,k,iOxyg),0.0_r8)
                fac3=fac2/(K_Nitri(ng)+fac2)
# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                fac1=fac*fac3*cff
# else
                fac1=fac*fac3
# endif
#else
                fac1=fac
#endif
                fac4=D_p5NH4(ng)+PAR1-2.0_r8*I_thNH4(ng)
                cff1=(PAR1-I_thNH4(ng))/MAX(fac4,eps)
                cff2=1.0_r8-MAX(0.0_r8,cff1)
                cff3=fac1*cff2
!               ========================================================
#ifdef OXYGEN
!>              Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
!>              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                    &
!>   &                            2.0_r8*tl_N_Flux_Nitrifi
                ad_N_Flux_Nitrifi=ad_N_Flux_Nitrifi-                    &
     &                            2.0_r8*ad_Bio(i,k,iOxyg)
#endif
!>              Bio2(i,k,iNO3_)=Bio(i,k,iNO3_)
!>              tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+tl_N_Flux_Nitrifi
                ad_N_Flux_Nitrifi=ad_N_Flux_Nitrifi+ad_Bio(i,k,iNO3_)

!>              tl_N_Flux_Nitrifi=tl_Bio(i,k,iNH4_)*cff3+               &
!>   &                            Bio(i,k,iNH4_)*tl_cff3
                ad_cff3=ad_cff3+Bio(i,k,iNH4_)*ad_N_Flux_Nitrifi
                ad_Bio(i,k,iNH4_)=ad_Bio(i,k,iNH4_)+                    &
     &                            ad_N_Flux_Nitrifi*cff3
                ad_N_Flux_Nitrifi=0.0_r8

!>              Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>              tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio(i,k,iNH4_)*    &
!>   &                             tl_cff3)/(1.0_r8+cff3)
                adfac=ad_Bio(i,k,iNH4_)/(1.0_r8+cff3)
                ad_cff3=ad_cff3-Bio(i,k,iNH4_)*adfac
                ad_Bio(i,k,iNH4_)=adfac

!>              tl_cff3=tl_fac1*cff2+fac1*tl_cff2
                ad_cff2=ad_cff2+fac1*ad_cff3
                ad_fac1=ad_fac1+ad_cff3*cff2
                ad_cff3=0.0_r8

!>              tl_cff2=-(0.5_r8-SIGN(0.5_r8,-cff1))*tl_cff1
                ad_cff1=ad_cff1-(0.5_r8-SIGN(0.5_r8,-cff1))*ad_cff2
                ad_cff2=0.0_r8

!>              tl_cff1=(tl_PAR-tl_I_thNH4-cff1*tl_fac4)/MAX(fac4,eps)
                adfac=ad_cff1/MAX(fac4,eps)
                ad_fac4=ad_fac4-cff1*adfac
                ad_I_thNH4=ad_I_thNH4-adfac
                ad_PAR=ad_PAR+adfac
                ad_cff1=0.0_r8

!>              tl_fac4=tl_D_p5NH4+tl_PAR-2.0_r8*tl_I_thNH4
                ad_I_thNH4=ad_I_thNH4-2.0_r8*ad_fac4
                ad_PAR=ad_PAR+ad_fac4
                ad_D_p5NH4=ad_D_p5NH4+ad_fac4
                ad_fac4=0.0_r8
#ifdef OXYGEN
# ifdef TDEPENDANCE
!>              tl_fac1=tl_fac*fac3*cff+                                &
!>   &                  fac*tl_fac3*cff+                                &
!>   &                  fac*fac3*tl_cff
                ad_cff=ad_cff+fac*fac3*ad_fac1
                ad_fac3=ad_fac3+fac*cff*ad_fac1
                ad_fac=ad_fac+fac3*cff*ad_fac1
                ad_fac1=0.0_r8

!>              tl_cff=cff*tl_Bio(i,k,itemp)*LOG(NitriR_t(ng))
#  ifndef UV_FIXED_TL
                ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+                    &
     &                            cff*ad_cff*LOG(NitriR_t(ng))
#  endif
                ad_cff=0.0_r8
# else
!>              tl_fac1=tl_fac*fac3+fac*tl_fac3
                ad_fac3=ad_fac3+fac*ad_fac1
                ad_fac=ad_fac+fac3*ad_fac1
                ad_fac1=0.0_r8
# endif
!>              tl_fac3=(tl_fac2-fac3*(tl_K_Nitri+tl_fac2))/            &
!>   &                  (K_Nitri(ng)+fac2)
                adfac=ad_fac3/(K_Nitri(ng)+fac2)
                ad_fac2=ad_fac2-fac3*adfac
                ad_K_Nitri=ad_K_Nitri-fac3*adfac
                ad_fac2=ad_fac2+adfac
                ad_fac3=0.0_r8

!>              tl_fac2=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)))*          &
!>   &                  tl_Bio(i,k,iOxyg)
                adfac=0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg))
                ad_Bio(i,k,iOxyg)=ad_Bio(i,k,iOxyg)+adfac*ad_fac2
                ad_fac2=0.0_r8
#else
!>              tl_fac1=tl_fac
                ad_fac=ad_fac+ad_fac1
                ad_fac1=0.0_r8
#endif
!>              tl_fac=dtdays*tl_NitriR
                ad_NitriR=ad_NitriR+dtdays*ad_fac
                ad_fac=0.0_r8
!=======================================================================
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                cff=PhyCN(ng)*12.0_r8

#ifdef PHYT2 /*1*/
                fac1=Bio1(i,k,iPhyt)*cff
                cff1=Bio1(i,k,iChlo1)/(fac1+eps)
                Chl2C1=MIN(cff1,Chl2C_m(ng))

                fac12=Bio1(i,k,iPhyt2)*cff
                cff12=Bio1(i,k,iChlo2)/(fac12+eps)
                Chl2C2=MIN(cff12,Chl2C_m(ng))
#else
                fac1=Bio1(i,k,iPhyt)*cff
                cff1=Bio1(i,k,iChlo)/(fac1+eps)
                Chl2C=MIN(cff1,Chl2C_m(ng))
#endif
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
# ifdef VP_TEMP1022
                Vp=Vp0(ng)*0.59_r8*(1.022_r8**Bio(i,k,itemp))
# elif defined VP_TEMP1033
                Vp=Vp0(ng)*0.59_r8*(1.033_r8**Bio(i,k,itemp))
# else
                Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
# endif
                fac1=PAR1*PhyIS(ng)
                fac=SQRT(Vp*Vp+fac1*fac1)
                Epp=Vp/fac
                t_PPmax=Epp*fac1
!
!  Nutrient-limitation terms -------------------------------------------
!
                cff1=Bio1(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio1(i,k,iNO3_)*K_NO3(ng)
                inhNH4=1.0_r8/(1.0_r8+cff1)
                L_NH4=cff1/(1.0_r8+cff1)
                fac2=cff2/(1.0_r8+cff2)
                L_NO3=fac2*inhNH4
                LTOT=L_NO3+L_NH4
#ifdef PHOSPHORUS
                cff3=Bio1(i,k,iPO4_)*K_PO4(ng)
                L_PO4=cff3/(1.0_r8+cff3)
                LMIN=MIN(LTOT,L_PO4)
#endif
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
                fac1=dtdays*t_PPmax*Bio1(i,k,iPhyt)
                fac4=fac1*K_NO3(ng)*inhNH4
                cff4=fac4/(1.0_r8+cff2)
                fac5=fac1*K_NH4(ng)
                cff5=fac5/(1.0_r8+cff1)
!!              Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
!!              Bio1(i,k,iNH4_)=Bio(i,k,iNH4_)
#ifdef PHOSPHORUS
                fac6=fac1*PhyPN(ng)*K_PO4(ng)
                cff6=fac6/(1.0_r8+cff3)
# ifdef PHYT2 /*1*/
                fac12=dtdays*t_PPmax2*Bio1(i,k,iPhyt2)
                cff42=fac12*K_NO3(ng)*inhNH4/(1.0_r8+cff2)
                cff52=fac12*K_NH4(ng)/(1.0_r8+cff1)
                fac62=fac12*PhyPN(ng)*K_PO4(ng)
                cff62=fac62/(1.0_r8+cff3)

# endif
!!              Bio1(i,k,iPO4_)=Bio(i,k,iPO4_)
                IF (LMIN.eq.L_PO4) THEN
!
!  Phosphorus-limitation for N_flux
!
!>                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff6)
                  P_Flux=Bio(i,k,iPO4_)*cff6
                  fac2=L_NO3/MAX(LTOT,eps)
                  fac=P_Flux/PhyPN(ng)
                  N_Flux_NewProd=fac*fac2
                  fac3=L_NH4/MAX(LTOT,eps)
                  N_Flux_RegProd=fac*fac3
# ifdef PHYT2 /*2*/

                  P_Flux2=Bio(i,k,iPO4_)*cff62
                  fac02=P_Flux2/PhyPN(ng)
                  N_Flux_NewProd2=fac02*fac2
                  N_Flux_RegProd2=fac02*fac3
# endif
!>                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd
!>                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd
                ELSE
!
!  Nitrogen-limitation for N_flux
!
!>                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
!>                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  N_Flux_NewProd=Bio2(i,k,iNO3_)*cff4
                  N_Flux_RegProd=Bio2(i,k,iNH4_)*cff5
                  P_Flux=(N_Flux_NewProd+N_Flux_RegProd)*PhyPN(ng)
# ifdef PHYT2 /*3*/

                  N_Flux_NewProd2=Bio2(i,k,iNO3_)*cff42
                  N_Flux_RegProd2=Bio2(i,k,iNH4_)*cff52
                  P_Flux2=(N_Flux_NewProd2+N_Flux_RegProd2)*PhyPN(ng)
# endif
!>                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-P_Flux
                END IF
#else
!
!  Nitrogen-limitation for N_flux
!
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                N_Flux_NewProd=Bio2(i,k,iNO3_)*cff4
                N_Flux_RegProd=Bio2(i,k,iNH4_)*cff5
#endif
!
!  Update phytoplankton ------------------------------------------------
!
!>              Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!>              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
!>   &                         N_Flux_NewProd+N_Flux_RegProd
!
!  Update chlorophyll
!
#ifdef PHOSPHORUS
                fac=dtdays*t_PPmax*t_PPmax*LMIN*LMIN
#else
                fac=dtdays*t_PPmax*t_PPmax*LTOT*LTOT
#endif

#ifdef PHYT2 /*4*/
                fac1=fac*Chl2C_m(ng)
                fac2=PhyIS(ng)*MAX(Chl2C1,eps)*PAR1
                fac3=fac1/(fac2+eps)

                fac02=dtdays*t_PPmax2*t_PPmax2*LMIN*LMIN
                fac12=fac02*Chl2C_m(ng)
                fac22=PhyIS(ng)*MAX(Chl2C2,eps)*PAR12
                fac32=fac12/(fac22+eps)
#else
                fac1=fac*Chl2C_m(ng)
                fac2=PhyIS(ng)*MAX(Chl2C,eps)*PAR1
                fac3=fac1/(fac2+eps)
#endif
!>              Bio1(i,k,iChlo)=Bio(i,k,iChlo)
!>              Bio(i,k,iChlo)=Bio(i,k,iChlo)+Bio(i,k,iChlo)*fac3
#ifdef OXYGEN
!
!  Update oxygen
!
!>              Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
!>              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
!>   &                         N_Flux_NewProd*rOxNO3+                   &
!>   &                         N_Flux_RegProd*rOxNH4
#endif
!               ========================================================
#ifdef OXYGEN
!
!  Update oxygen
!
#ifdef PHYT2 /*5*/
!>              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)                     &
!>   &                            +tl_N_Flux_NewProd2*rOxNO3            &
!>   &                            +tl_N_Flux_RegProd2*rOxNH4            &
!>   &                            +tl_N_Flux_NewProd*rOxNO3             &
!>   &                            +tl_N_Flux_RegProd*rOxNH4
                adfac=ad_Bio(i,k,iOxyg)
                ad_N_Flux_RegProd2=ad_N_Flux_RegProd2+adfac*rOxNH4
                ad_N_Flux_NewProd2=ad_N_Flux_NewProd2+adfac*rOxNO3

#endif
!>              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)+                    &
!>   &                            tl_N_Flux_NewProd*rOxNO3+             &
!>   &                            tl_N_Flux_RegProd*rOxNH4
                adfac=ad_Bio(i,k,iOxyg)
                ad_N_Flux_RegProd=ad_N_Flux_RegProd+adfac*rOxNH4
                ad_N_Flux_NewProd=ad_N_Flux_NewProd+adfac*rOxNO3
#endif
!
!  Update chlorophyll
!
#ifdef PHYT2 /*6*/
!>              tl_Bio(i,k,iChlo2)=tl_Bio(i,k,iChlo2)+                  &
!>   &                             tl_Bio(i,k,iChlo2)*fac32+            &
!>   &                             Bio1(i,k,iChlo2)*tl_fac32
                ad_fac32=ad_fac32+Bio1(i,k,iChlo2)*ad_Bio(i,k,iChlo2)
                ad_Bio(i,k,iChlo2)=ad_Bio(i,k,iChlo2)*(1.0_r8+fac12) 

!>              tl_Bio(i,k,iChlo1)=tl_Bio(i,k,iChlo1)+                  &
!>   &                            tl_Bio(i,k,iChlo1)*fac3+              &
!>   &                            Bio1(i,k,iChlo1)*tl_fac3
                ad_fac3=ad_fac3+Bio1(i,k,iChlo1)*ad_Bio(i,k,iChlo1)
                ad_Bio(i,k,iChlo1)=ad_Bio(i,k,iChlo1)*(1.0_r8+fac1) 

!>              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)+                    &
!>   &                            tl_Bio(i,k,iChlo1)*fac3+              &
!>   &                            Bio1(i,k,iChlo1)*tl_fac3+             &
!>   &                            tl_Bio(i,k,iChlo2)*fac32+             &
!>   &                            Bio1(i,k,iChlo2)*tl_fac32
!                ad_fac3=ad_fac3+Bio1(i,k,iChlo1)*ad_Bio(i,k,iChlo)
!                ad_fac32=ad_fac32+Bio1(i,k,iChlo2)*ad_Bio(i,k,iChlo)
!                ad_Bio(i,k,iChlo1)=ad_Bio(i,k,iChlo1)+fac1*ad_Bio(i,k,iChlo)
!                ad_Bio(i,k,iChlo2)=ad_Bio(i,k,iChlo2)+fac12*ad_Bio(i,k,iChlo)
!                ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)+0.5*(fac1+fac12)*ad_Bio(i,k,iChlo)


!>              tl_fac32=(tl_fac12-tl_fac22*fac32)/(fac22+eps)
                adfac=ad_fac32/(fac22+eps)
                ad_fac22=ad_fac22-fac32*adfac
                ad_fac12=ad_fac12+adfac
                ad_fac32=0.0_r8

!>              tl_fac22=tl_PhyIS*MAX(Chl2C2,eps)*PAR2+                 &
!>   &                   PhyIS(ng)*tl_Chl2C2*PAR2+                      &
!>   &                   PhyIS(ng)*MAX(Chl2C2,eps)*tl_PAR2
                ad_PAR2=ad_PAR2+PhyIS(ng)*MAX(Chl2C2,eps)*ad_fac22
                ad_Chl2C2=ad_Chl2C2+PhyIS(ng)*PAR12*ad_fac22
                ad_PhyIS=ad_PhyIS+MAX(Chl2C2,eps)*PAR12*ad_fac22
                ad_fac22=0.0_r8

!>              tl_fac12=tl_fac02*Chl2C_m(ng)+fac02*tl_Chl2C_m
                ad_fac02=ad_fac02+Chl2C_m(ng)*ad_fac12
                ad_Chl2C_m=ad_Chl2C_m+fac02*ad_fac12
                ad_fac12=0.0_r8

!>              tl_fac02=dtdays*2.0_r8*(t_PPmax2*tl_t_PPmax2*LMIN*LMIN+ &
!>   &                                  t_PPmax2*t_PPmax2*LMIN*tl_LMIN)
                adfac=dtdays*2.0_r8*t_PPmax2*LMIN*ad_fac02
                ad_LMIN=ad_LMIN+t_PPmax2*adfac
                ad_t_PPmax2=ad_t_PPmax2+LMIN*adfac
                ad_fac=0.0_r8

!>              tl_fac3=(tl_fac1-tl_fac2*fac3)/(fac2+eps)
                adfac=ad_fac3/(fac2+eps)
                ad_fac2=ad_fac2-fac3*adfac
                ad_fac1=ad_fac1+adfac
                ad_fac3=0.0_r8

!>              tl_fac2=tl_PhyIS*MAX(Chl2C,eps)*PAR+                    &
!>   &                  PhyIS(ng)*tl_Chl2C*PAR+                         &
!>   &                  PhyIS(ng)*MAX(Chl2C,eps)*tl_PAR
                ad_PAR=ad_PAR+PhyIS(ng)*MAX(Chl2C1,eps)*ad_fac2
                ad_Chl2C1=ad_Chl2C1+PhyIS(ng)*PAR1*ad_fac2
                ad_PhyIS=ad_PhyIS+MAX(Chl2C1,eps)*PAR1*ad_fac2
                ad_fac2=0.0_r8

!>              tl_fac1=tl_fac*Chl2C_m(ng)+fac*tl_Chl2C_m
                ad_fac=ad_fac+Chl2C_m(ng)*ad_fac1
                ad_Chl2C_m=ad_Chl2C_m+fac*ad_fac1
                ad_fac1=0.0_r8
#else /*6*/
!>              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)+                    &
!>   &                            tl_Bio(i,k,iChlo)*fac3+               &
!>   &                            Bio1(i,k,iChlo)*tl_fac3
                ad_fac3=ad_fac3+Bio1(i,k,iChlo)*ad_Bio(i,k,iChlo)
                ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)*(1.0_r8+fac1) /**Noda**/

!>              tl_fac3=(tl_fac1-tl_fac2*fac3)/(fac2+eps)
                adfac=ad_fac3/(fac2+eps)
                ad_fac2=ad_fac2-fac3*adfac
                ad_fac1=ad_fac1+adfac
                ad_fac3=0.0_r8

!>              tl_fac2=tl_PhyIS*MAX(Chl2C,eps)*PAR+                    &
!>   &                  PhyIS(ng)*tl_Chl2C*PAR+                         &
!>   &                  PhyIS(ng)*MAX(Chl2C,eps)*tl_PAR
                ad_PAR=ad_PAR+PhyIS(ng)*MAX(Chl2C,eps)*ad_fac2
                ad_Chl2C=ad_Chl2C+PhyIS(ng)*PAR1*ad_fac2
                ad_PhyIS=ad_PhyIS+MAX(Chl2C,eps)*PAR1*ad_fac2
                ad_fac2=0.0_r8

!>              tl_fac1=tl_fac*Chl2C_m(ng)+fac*tl_Chl2C_m
                ad_fac=ad_fac+Chl2C_m(ng)*ad_fac1
                ad_Chl2C_m=ad_Chl2C_m+fac*ad_fac1
                ad_fac1=0.0_r8
#endif /*6*/
#if defined PHOSPHORUS
!>              tl_fac=dtdays*2.0_r8*(t_PPmax*tl_t_PPmax*LMIN*LMIN+     &
!>   &                                t_PPmax*t_PPmax*LMIN*tl_LMIN)
                adfac=dtdays*2.0_r8*t_PPmax*LMIN*ad_fac
                ad_LMIN=ad_LMIN+t_PPmax*adfac
                ad_t_PPmax=ad_t_PPmax+LMIN*adfac
                ad_fac=0.0_r8
#else
!>              tl_fac=dtdays*2.0_r8*(t_PPmax*tl_t_PPmax*LTOT*LTOT+     &
!>   &                                t_PPmax*t_PPmax*LTOT*tl_LTOT)
                adfac=dtdays*2.0_r8*t_PPmax*LTOT*ad_fac
                ad_LTOT=ad_LTOT+t_PPmax*adfac
                ad_t_PPmax=ad_t_PPmax+LTOT*adfac
                ad_fac=0.0_r8
#endif
!
!  Update phytoplankton
!

#ifdef PHYT2 /*7*/
!>              tl_Bio(i,k,iPhyt2)=tl_Bio(i,k,iPhyt2)+                  &
!>   &                            tl_N_Flux_NewProd2+tl_N_Flux_RegProd2
                ad_N_Flux_RegProd2=ad_N_Flux_RegProd2+ad_Bio(i,k,iPhyt2)
                ad_N_Flux_NewProd2=ad_N_Flux_NewProd2+ad_Bio(i,k,iPhyt2)

#endif
!>              Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!>              tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+                    &
!>   &                            tl_N_Flux_NewProd+tl_N_Flux_RegProd
                ad_N_Flux_RegProd=ad_N_Flux_RegProd+ad_Bio(i,k,iPhyt)
                ad_N_Flux_NewProd=ad_N_Flux_NewProd+ad_Bio(i,k,iPhyt)
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
!               ========================================================
                cff1=Bio1(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio1(i,k,iNO3_)*K_NO3(ng)
                fac2=cff2/(1.0_r8+cff2)
#ifdef PHOSPHORUS
                cff3=Bio1(i,k,iPO4_)*K_PO4(ng)
#endif
                fac1=dtdays*t_PPmax*Bio1(i,k,iPhyt)
                fac4=fac1*K_NO3(ng)*inhNH4
                cff4=fac4/(1.0_r8+cff2)
                fac5=fac1*K_NH4(ng)
                cff5=fac5/(1.0_r8+cff1)
#ifdef PHOSPHORUS
                fac6=fac1*PhyPN(ng)*K_PO4(ng)
                cff6=fac6/(1.0_r8+cff3)

# ifdef PHYT2 /*8*/
                fac12=dtdays*t_PPmax2*Bio1(i,k,iPhyt2)
                fac42=fac12*K_NO3(ng)*inhNH4
                cff42=fac42/(1.0_r8+cff2)
                fac52=fac12*K_NH4(ng)
                cff52=fac52/(1.0_r8+cff1)
                fac62=fac12*PhyPN(ng)*K_PO4(ng)
                cff62=fac62/(1.0_r8+cff3)
# endif
                IF (LMIN.eq.L_PO4) THEN
                  fac2=L_NO3/MAX(LTOT,eps)
                  fac=P_Flux/PhyPN(ng)
                  fac3=L_NH4/MAX(LTOT,eps)

# ifdef PHYT2 /*9*/
                  fac02=P_Flux2/PhyPN(ng)
# endif
                END IF
#endif
!               ========================================================
#ifdef PHOSPHORUS
                IF (LMIN.eq.L_PO4) THEN
!
!  Phosphorus-limitation for N_flux
!

# ifdef PHYT2 /*10*/
!>                tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)                   &
!>   &                              -tl_N_Flux_RegProd2                 &
!>   &                              -tl_N_Flux_RegProd
                  ad_N_Flux_RegProd2=ad_N_Flux_RegProd2-ad_Bio(i,k,iNH4_)
                  ad_N_Flux_RegProd=ad_N_Flux_RegProd-ad_Bio(i,k,iNH4_)

# else
!>                tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)-tl_N_Flux_RegProd
                  ad_N_Flux_RegProd=ad_N_Flux_RegProd-ad_Bio(i,k,iNH4_)

# endif

# ifdef PHYT2 /*11*/
!>                tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)                   &
!>   &                              -tl_N_Flux_NewProd2                 &
!>   &                              -tl_N_Flux_NewProd
                  ad_N_Flux_NewProd2=ad_N_Flux_NewProd2-ad_Bio(i,k,iNO3_)
                  ad_N_Flux_NewProd=ad_N_Flux_NewProd-ad_Bio(i,k,iNO3_)

!>                tl_N_Flux_RegProd2=tl_fac02*fac3+fac02*tl_fac3
                  ad_fac3=ad_fac3+fac02*ad_N_Flux_RegProd2
                  ad_fac02=ad_fac02+fac3*ad_N_Flux_RegProd2
                  ad_N_Flux_RegProd2=0.0_r8

!>                tl_N_Flux_NewProd2=tl_fac02*fac2+fac02*tl_fac2
                  ad_fac2=ad_fac2+fac02*ad_N_Flux_NewProd2
                  ad_fac02=ad_fac02+fac2*ad_N_Flux_NewProd2
                  ad_N_Flux_NewProd2=0.0_r8

!>                tl_fac02=(tl_P_Flux2-fac02*tl_PhyPN)/PhyPN(ng)
                  adfac=ad_fac02/PhyPN(ng)
                  ad_PhyPN=ad_PhyPN-fac02*adfac
                  ad_P_Flux2=ad_P_Flux2+adfac
                  ad_fac02=0.0_r8

!>                tl_P_Flux2=tl_Bio(i,k,iPO4_)*cff62+                   &
!>   &                      Bio(i,k,iPO4_)*tl_cff62
                  ad_cff62=ad_cff62+Bio(i,k,iPO4_)*ad_P_Flux2
                  ad_Bio(i,k,iPO4_)=ad_Bio(i,k,iPO4_)+ad_P_Flux2*cff62
                  ad_P_Flux2=0.0_r8
# else 
!>                tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)-tl_N_Flux_NewProd
                  ad_N_Flux_NewProd=ad_N_Flux_NewProd-ad_Bio(i,k,iNO3_)
# endif

!>                tl_N_Flux_RegProd=tl_fac*fac3+fac*tl_fac3
                  ad_fac3=ad_fac3+fac*ad_N_Flux_RegProd
                  ad_fac=ad_fac+fac3*ad_N_Flux_RegProd
                  ad_N_Flux_RegProd=0.0_r8

!>                tl_fac3=(tl_L_NH4-fac3*tl_LTOT)/MAX(LTOT,eps)
                  adfac=ad_fac3/MAX(LTOT,eps)
                  ad_LTOT=ad_LTOT-fac3*adfac
                  ad_L_NH4=ad_L_NH4+adfac
                  ad_fac3=0.0_r8

!>                tl_N_Flux_NewProd=tl_fac*fac2+fac*tl_fac2
                  ad_fac2=ad_fac2+fac*ad_N_Flux_NewProd
                  ad_fac=ad_fac+fac2*ad_N_Flux_NewProd
                  ad_N_Flux_NewProd=0.0_r8

!>                tl_fac=(tl_P_Flux-fac*tl_PhyPN)/PhyPN(ng)
                  adfac=ad_fac/PhyPN(ng)
                  ad_PhyPN=ad_PhyPN-fac*adfac
                  ad_P_Flux=ad_P_Flux+adfac
                  ad_fac=0.0_r8

!>                tl_fac2=(tl_L_NO3-fac2*tl_LTOT)/MAX(LTOT,eps)
                  adfac=ad_fac2/MAX(LTOT,eps)
                  ad_LTOT=ad_LTOT-fac2*adfac
                  ad_L_NO3=ad_L_NO3+adfac
                  ad_fac2=0.0_r8

!>                tl_P_Flux=tl_Bio(i,k,iPO4_)*cff6+                     &
!>   &                      Bio(i,k,iPO4_)*tl_cff6
                  ad_cff6=ad_cff6+Bio(i,k,iPO4_)*ad_P_Flux
                  ad_Bio(i,k,iPO4_)=ad_Bio(i,k,iPO4_)+ad_P_Flux*cff6
                  ad_P_Flux=0.0_r8

# ifdef PHYT2 /*12*/
!>                tl_Bio(i,k,iPO4_)=(tl_Bio(i,k,iPO4_)-Bio(i,k,iPO4_)*  &
!>   &                               (tl_cff6+tl_cff62))/(1.0_r8+cff_)
                  cff_=cff6+cff62
                  adfac=ad_Bio(i,k,iPO4_)/(1.0_r8+cff_)
                  ad_cff6=ad_cff6-Bio(i,k,iPO4_)*adfac
                  ad_cff62=ad_cff62-Bio(i,k,iPO4_)*adfac
                  ad_Bio(i,k,iPO4_)=adfac
# else
!>                tl_Bio(i,k,iPO4_)=(tl_Bio(i,k,iPO4_)-Bio(i,k,iPO4_)*  &
!>   &                               tl_cff6)/(1.0_r8+cff6)
                  adfac=ad_Bio(i,k,iPO4_)/(1.0_r8+cff6)
                  ad_cff6=ad_cff6-Bio(i,k,iPO4_)*adfac
                  ad_Bio(i,k,iPO4_)=adfac
# endif
                ELSE
!
!  Nitrogen-limitation for N_flux
!

# ifdef PHYT2 /*13*/
!>                tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)                   &
!>   &                              -tl_P_Flux2                         &
!>   &                              -tl_P_Flux
                  ad_P_Flux2=ad_P_Flux2-ad_Bio(i,k,iPO4_)
                  ad_P_Flux=ad_P_Flux-ad_Bio(i,k,iPO4_)

!>                tl_P_Flux2=(tl_N_Flux_NewProd2+tl_N_Flux_RegProd2)*   &
!>   &                      PhyPN(ng)+                                  &
!>   &                      (N_Flux_NewProd2+N_Flux_RegProd2)*tl_PhyPN
                  adfac=ad_P_Flux2*PhyPN(ng)
                  ad_N_Flux_RegProd2=ad_N_Flux_RegProd2+adfac
                  ad_N_Flux_NewProd2=ad_N_Flux_NewProd2+adfac
                  ad_PhyPN=ad_PhyPN+(N_Flux_NewProd2+N_Flux_RegProd2)*  &
     &                              ad_P_Flux2
                  ad_P_Flux2=0.0_r8

!>                tl_N_Flux_RegProd2=tl_Bio(i,k,iNH4_)*cff52+           &
!>   &                              Bio2(i,k,iNH4_)*tl_cff52
                  adfac=ad_N_Flux_RegProd2
                  ad_cff52=ad_cff52+Bio2(i,k,iNH4_)*adfac
                  ad_Bio(i,k,iNH4_)=ad_Bio(i,k,iNH4_)+adfac*cff52
                  ad_N_Flux_RegProd2=0.0_r8

!>                tl_N_Flux_NewProd2=tl_Bio(i,k,iNO3_)*cff42+           &
!>   &                              Bio2(i,k,iNO3_)*tl_cff42
                  adfac=ad_N_Flux_NewProd2
                  ad_cff42=ad_cff42+Bio2(i,k,iNO3_)*adfac
                  ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)+adfac*cff42
                  ad_N_Flux_NewProd2=0.0_r8

# else
!>                tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)-tl_P_Flux
                  ad_P_Flux=ad_P_Flux-ad_Bio(i,k,iPO4_)

# endif
!>                tl_P_Flux=(tl_N_Flux_NewProd+tl_N_Flux_RegProd)*      &
!>   &                      PhyPN(ng)+                                  &
!>   &                      (N_Flux_NewProd+N_Flux_RegProd)*tl_PhyPN
                  adfac=ad_P_Flux*PhyPN(ng)
                  ad_N_Flux_RegProd=ad_N_Flux_RegProd+adfac
                  ad_N_Flux_NewProd=ad_N_Flux_NewProd+adfac
                  ad_PhyPN=ad_PhyPN+(N_Flux_NewProd+N_Flux_RegProd)*    &
     &                              ad_P_Flux
                  ad_P_Flux=0.0_r8

!>                tl_N_Flux_RegProd=tl_Bio(i,k,iNH4_)*cff5+             &
!>   &                              Bio2(i,k,iNH4_)*tl_cff5
                  adfac=ad_N_Flux_RegProd
                  ad_cff5=ad_cff5+Bio2(i,k,iNH4_)*adfac
                  ad_Bio(i,k,iNH4_)=ad_Bio(i,k,iNH4_)+adfac*cff5
                  ad_N_Flux_RegProd=0.0_r8

!>                tl_N_Flux_NewProd=tl_Bio(i,k,iNO3_)*cff4+             &
!>   &                              Bio2(i,k,iNO3_)*tl_cff4
                  adfac=ad_N_Flux_NewProd
                  ad_cff4=ad_cff4+Bio2(i,k,iNO3_)*adfac
                  ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)+adfac*cff4
                  ad_N_Flux_NewProd=0.0_r8

# ifdef PHYT2 /*14*/
!>                tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio2(i,k,iNH4_)* &
!>   &                               (tl_cff5+tl_cff52))/(1.0_r8+cff_)
                  cff_=cff5+cff52
                  adfac=ad_Bio(i,k,iNH4_)/(1.0_r8+cff_)
                  ad_cff52=ad_cff52-Bio2(i,k,iNH4_)*adfac
                  ad_cff5=ad_cff5-Bio2(i,k,iNH4_)*adfac
                  ad_Bio(i,k,iNH4_)=adfac

!>                tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio2(i,k,iNO3_)* &
!>   &                               (tl_cff4+tl_cff42))/(1.0_r8+cff_)
                  cff_=cff4+cff42
                  adfac=ad_Bio(i,k,iNO3_)/(1.0_r8+cff_)
                  ad_cff42=ad_cff42-Bio2(i,k,iNO3_)*adfac
                  ad_cff4=ad_cff4-Bio2(i,k,iNO3_)*adfac
                  ad_Bio(i,k,iNO3_)=adfac

# else
!>                tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio2(i,k,iNH4_)* &
!>   &                               tl_cff5)/(1.0_r8+cff5)
                  adfac=ad_Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  ad_cff5=ad_cff5-Bio2(i,k,iNH4_)*adfac
                  ad_Bio(i,k,iNH4_)=adfac

!>                tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio2(i,k,iNO3_)* &
!>   &                               tl_cff4)/(1.0_r8+cff4)
                  adfac=ad_Bio(i,k,iNO3_)/(1.0_r8+cff4)
                  ad_cff4=ad_cff4-Bio2(i,k,iNO3_)*adfac
                  ad_Bio(i,k,iNO3_)=adfac
# endif
                END IF

# ifdef PHYT2 /*15*/
!>              tl_cff62=(tl_fac62-cff62*tl_cff3)/(1.0_r8+cff3)
                adfac=ad_cff62/(1.0_r8+cff3)
                ad_cff3=ad_cff3-cff62*adfac
                ad_fac62=ad_fac62+adfac
                ad_cff62=0.0_r8

!>              tl_fac62=tl_fac12*PhyPN(ng)*K_PO4(ng)+                  &
!>   &                  fac12*tl_PhyPN*K_PO4(ng)+                       &
!>   &                  fac12*PhyPN(ng)*tl_K_PO4
                ad_K_PO4=ad_K_PO4+fac12*PhyPN(ng)*ad_fac62
                ad_PhyPN=ad_PhyPN+fac12*K_PO4(ng)*ad_fac62
                ad_fac12=ad_fac12+PhyPN(ng)*K_PO4(ng)*ad_fac62
                ad_fac62=0.0_r8

!>              tl_cff52=(tl_fac52-cff52*tl_cff1)/(1.0_r8+cff1)
                adfac=ad_cff52/(1.0_r8+cff1)
                ad_cff1=ad_cff1-cff52*adfac
                ad_fac52=ad_fac52+adfac
                ad_cff52=0.0_r8

!>              tl_fac52=tl_fac12*K_NH4(ng)+fac12*tl_K_NH4
                ad_K_NH4=ad_K_NH4+fac12*ad_fac52
                ad_fac12=ad_fac12+K_NH4(ng)*ad_fac52
                ad_fac52=0.0_r8

!>              tl_cff42=(tl_fac42-cff42*tl_cff2)/(1.0_r8+cff2)
                adfac=ad_cff42/(1.0_r8+cff2)
                ad_cff2=ad_cff2-cff42*adfac
                ad_fac42=ad_fac42+adfac
                ad_cff42=0.0_r8

!>              tl_fac42=tl_fac12*K_NO3(ng)*inhNH4+                     &
!>   &                  fac12*tl_K_NO3*inhNH4+                          &
!>   &                  fac12*K_NO3(ng)*tl_inhNH4
                ad_inhNH4=ad_inhNH4+fac12*K_NO3(ng)*ad_fac42
                ad_K_NO3=ad_K_NO3+fac12*inhNH4*ad_fac42
                ad_fac12=ad_fac12+K_NO3(ng)*inhNH4*ad_fac42
                ad_fac42=0.0_r8

!>              tl_fac12=dtdays*(tl_t_PPmax2*Bio1(i,k,iPhyt2)+          &
!>   &                          t_PPmax2*tl_Bio(i,k,iPhyt2))
                adfac=dtdays*ad_fac12
                ad_Bio(i,k,iPhyt2)=ad_Bio(i,k,iPhyt2)+t_PPmax2*adfac
                ad_t_PPmax2=ad_t_PPmax2+Bio1(i,k,iPhyt2)*adfac
                ad_fac12=0.0_r8

# endif 
!>              tl_cff6=(tl_fac6-cff6*tl_cff3)/(1.0_r8+cff3)
                adfac=ad_cff6/(1.0_r8+cff3)
                ad_cff3=ad_cff3-cff6*adfac
                ad_fac6=ad_fac6+adfac
                ad_cff6=0.0_r8

!>              tl_fac6=tl_fac1*PhyPN(ng)*K_PO4(ng)+                    &
!>   &                  fac1*tl_PhyPN*K_PO4(ng)+                        &
!>   &                  fac1*PhyPN(ng)*tl_K_PO4
                ad_K_PO4=ad_K_PO4+fac1*PhyPN(ng)*ad_fac6
                ad_PhyPN=ad_PhyPN+fac1*K_PO4(ng)*ad_fac6
                ad_fac1=ad_fac1+PhyPN(ng)*K_PO4(ng)*ad_fac6
                ad_fac6=0.0_r8
#else
!
!  Nitrogen-limitation for N_flux
!
!>              tl_N_Flux_RegProd=tl_Bio(i,k,iNH4_)*cff5+               &
!>   &                            Bio2(i,k,iNH4_)*tl_cff5
                adfac=ad_N_Flux_RegProd
                ad_cff5=ad_cff5+Bio2(i,k,iNH4_)*adfac
                ad_Bio(i,k,iNH4_)=ad_Bio(i,k,iNH4_)+adfac*cff5
                ad_N_Flux_RegProd=0.0_r8

!>              tl_N_Flux_NewProd=tl_Bio(i,k,iNO3_)*cff4+               &
!>   &                            Bio2(i,k,iNO3_)*tl_cff4
                adfac=ad_N_Flux_NewProd
                ad_cff4=ad_cff4+Bio2(i,k,iNO3_)*adfac
                ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)+adfac*cff4
                ad_N_Flux_NewProd=0.0_r8

!>              tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio2(i,k,iNH4_)*   &
!>   &                             tl_cff5)/(1.0_r8+cff5)
                adfac=ad_Bio(i,k,iNH4_)/(1.0_r8+cff5)
                ad_cff5=ad_cff5-Bio2(i,k,iNH4_)*adfac
                ad_Bio(i,k,iNH4_)=adfac

!>              tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio2(i,k,iNO3_)*   &
!>   &                             tl_cff4)/(1.0_r8+cff4)
                adfac=ad_Bio(i,k,iNO3_)/(1.0_r8+cff4)
                ad_cff4=ad_cff4-Bio2(i,k,iNO3_)*adfac
                ad_Bio(i,k,iNO3_)=adfac
#endif
!>              tl_cff5=(tl_fac5-cff5*tl_cff1)/(1.0_r8+cff1)
                adfac=ad_cff5/(1.0_r8+cff1)
                ad_cff1=ad_cff1-cff5*adfac
                ad_fac5=ad_fac5+adfac
                ad_cff5=0.0_r8

!>              tl_fac5=tl_fac1*K_NH4(ng)+fac1*tl_K_NH4
                ad_K_NH4=ad_K_NH4+fac1*ad_fac5
                ad_fac1=ad_fac1+K_NH4(ng)*ad_fac5
                ad_fac5=0.0_r8

!>              tl_cff4=(tl_fac4-cff4*tl_cff2)/(1.0_r8+cff2)
                adfac=ad_cff4/(1.0_r8+cff2)
                ad_cff2=ad_cff2-cff4*adfac
                ad_fac4=ad_fac4+adfac
                ad_cff4=0.0_r8

!>              tl_fac4=tl_fac1*K_NO3(ng)*inhNH4+                       &
!>   &                  fac1*tl_K_NO3*inhNH4+                           &
!>   &                  fac1*K_NO3(ng)*tl_inhNH4
                ad_inhNH4=ad_inhNH4+fac1*K_NO3(ng)*ad_fac4
                ad_K_NO3=ad_K_NO3+fac1*inhNH4*ad_fac4
                ad_fac1=ad_fac1+K_NO3(ng)*inhNH4*ad_fac4
                ad_fac4=0.0_r8

!>              tl_fac1=dtdays*(tl_t_PPmax*Bio1(i,k,iPhyt)+             &
!>   &                          t_PPmax*tl_Bio(i,k,iPhyt))
                adfac=dtdays*ad_fac1
                ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+t_PPmax*adfac
                ad_t_PPmax=ad_t_PPmax+Bio1(i,k,iPhyt)*adfac
                ad_fac1=0.0_r8
!
!  Nutrient-limitation terms
!
!               ========================================================
                cff1=Bio1(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio1(i,k,iNO3_)*K_NO3(ng)
                fac2=cff2/(1.0_r8+cff2)
#ifdef PHOSPHORUS
                cff3=Bio1(i,k,iPO4_)*K_PO4(ng)
#endif
!               ========================================================
#ifdef PHOSPHORUS
!>              tlfac=SIGN(0.5_r8,L_PO4-LTOT)
!>              tl_LMIN=(0.5_r8+tlfac)*tl_LTOT+(0.5_r8-tlfac)*tl_L_PO4
                adfac=SIGN(0.5_r8,L_PO4-LTOT)
                ad_L_PO4=ad_L_PO4+(0.5_r8-adfac)*ad_LMIN
                ad_LTOT=ad_LTOT+(0.5_r8+adfac)*ad_LMIN
                ad_LMIN=0.0_r8

!>              tl_L_PO4=tl_cff3*(1.0_r8-L_PO4)/(1.0_r8+cff3)
                ad_cff3=ad_cff3+ad_L_PO4*(1.0_r8-L_PO4)/(1.0_r8+cff3)
                ad_L_PO4=0.0_r8

!>              tl_cff3=tl_Bio(i,k,iPO4_)*K_PO4(ng)+                    &
!>   &                  Bio1(i,k,iPO4_)*tl_K_PO4
                ad_Bio(i,k,iPO4_)=ad_Bio(i,k,iPO4_)+ad_cff3*K_PO4(ng)
                ad_K_PO4=ad_K_PO4+Bio1(i,k,iPO4_)*ad_cff3
                ad_cff3=0.0_r8
#endif
!>              tl_LTOT=tl_L_NO3+tl_L_NH4
                ad_L_NH4=ad_L_NH4+ad_LTOT
                ad_L_NO3=ad_L_NO3+ad_LTOT
                ad_LTOT=0.0_r8

!>              tl_L_NO3=tl_fac2*inhNH4+fac2*tl_inhNH4
                ad_inhNH4=ad_inhNH4+fac2*ad_L_NO3
                ad_fac2=ad_fac2+inhNH4*ad_L_NO3
                ad_L_NO3=0.0_r8

!>              tl_fac2=tl_cff2*(1.0_r8-fac2)/(1.0_r8+cff2)
                ad_cff2=ad_cff2+ad_fac2*(1.0_r8-fac2)/(1.0_r8+cff2)
                ad_fac2=0.0_r8

!>              tl_L_NH4=tl_cff1*(1.0_r8-L_NH4)/(1.0_r8+cff1)
                ad_cff1=ad_cff1+ad_L_NH4*(1.0_r8-L_NH4)/(1.0_r8+cff1)
                ad_L_NH4=0.0_r8

!>              tl_inhNH4=-inhNH4*inhNH4*tl_cff1
                ad_cff1=ad_cff1-inhNH4*inhNH4*ad_inhNH4
                ad_inhNH4=0.0_r8

!>              tl_cff2=tl_Bio(i,k,iNO3_)*K_NO3(ng)+                    &
!>   &                  Bio1(i,k,iNO3_)*tl_K_NO3
                ad_Bio(i,k,iNO3_)=ad_Bio(i,k,iNO3_)+ad_cff2*K_NO3(ng)
                ad_K_NO3=ad_K_NO3+Bio1(i,k,iNO3_)*ad_cff2
                ad_cff2=0.0_r8

!>              tl_cff1=tl_Bio(i,k,iNH4_)*K_NH4(ng)+                    &
!>   &                  Bio1(i,k,iNH4_)*tl_K_NH4
                ad_Bio(i,k,iNH4_)=ad_Bio(i,k,iNH4_)+ad_cff1*K_NH4(ng)
                ad_K_NH4=ad_K_NH4+Bio1(i,k,iNH4_)*ad_cff1
                ad_cff1=0.0_r8


!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
!               ========================================================
                fac1=PAR1*PhyIS(ng)
                fac=SQRT(Vp*Vp+fac1*fac1)
!               ========================================================

!>              tl_t_PPmax=tl_Epp*fac1+Epp*tl_fac1
                ad_fac1=ad_fac1+Epp*ad_t_PPmax
                ad_Epp=ad_Epp+fac1*ad_t_PPmax
                ad_t_PPmax=0.0_r8

!>              tl_Epp=(tl_Vp-Epp*tl_fac)/fac
                adfac=ad_Epp/fac
                ad_fac=ad_fac-Epp*adfac
                ad_Vp=ad_Vp+adfac
                ad_Epp=0.0_r8

!>              tl_fac=(Vp*tl_Vp+fac1*tl_fac1)/fac
                adfac=ad_fac/fac
                ad_fac1=ad_fac1+fac1*adfac
                ad_Vp=ad_Vp+Vp*adfac
                ad_fac=0.0_r8

!>              tl_fac1=tl_PAR*PhyIS(ng)+PAR*tl_PhyIS
                ad_PAR=ad_PAR+ad_fac1*PhyIS(ng)
                ad_PhyIS=ad_PhyIS+PAR1*ad_fac1
                ad_fac1=0.0_r8

#ifdef GROWTH1 /*Noda*/

#else /*Noda*/
# ifdef VP_TEMP1022
!>              tl_Vp=tl_Vp0*0.59_r8*(1.022_r8**Bio(i,k,itemp))+        &
!>   &                Vp*tl_Bio(i,k,itemp)
#  ifndef UV_FIXED_TL
                ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+Vp*ad_Vp
#  endif
                ad_Vp0=ad_Vp0+ad_Vp*0.59_r8*(1.022_r8**Bio(i,k,itemp))
                ad_Vp=0.0_r8
# elif defined VP_TEMP1033
!>              tl_Vp=tl_Vp0*0.59_r8*(1.033_r8**Bio(i,k,itemp))+        &
!>   &                Vp*tl_Bio(i,k,itemp)
#  ifndef UV_FIXED_TL
                ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+Vp*ad_Vp
#  endif
                ad_Vp0=ad_Vp0+ad_Vp*0.59_r8*(1.033_r8**Bio(i,k,itemp))
                ad_Vp=0.0_r8
# else
!>              tl_Vp=tl_Vp0*0.59_r8*(1.066_r8**Bio(i,k,itemp))+        &
!>   &                Vp*tl_Bio(i,k,itemp)
#  ifndef UV_FIXED_TL
                ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+Vp*ad_Vp
#  endif
                ad_Vp0=ad_Vp0+ad_Vp*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                ad_Vp=0.0_r8
# endif
#endif /*Noda*/
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
!               ========================================================
                cff=PhyCN(ng)*12.0_r8
                fac1=Bio1(i,k,iPhyt)*cff

# ifdef PHYT2 /*1*/
                cff1=Bio1(i,k,iChlo1)/(fac1+eps)
                fac12=Bio1(i,k,iPhyt2)*cff
                cff12=Bio1(i,k,iChlo2)/(fac12+eps)
# else
                cff1=Bio1(i,k,iChlo)/(fac1+eps)
# endif
!               ========================================================

# ifdef PHYT2 /*2*/
!>              tlfac02=SIGN(0.5_r8,Chl2C_m(ng)-cff12)
!>              tl_Chl2C2=(0.5_r8+tlfac02)*tl_cff12+                    &
!>   &                    (0.5_r8-tlfac02)*tl_Chl2C_m
                adfac=SIGN(0.5_r8,Chl2C_m(ng)-cff12)
                ad_Chl2C_m=ad_Chl2C_m+(0.5_r8-adfac)*ad_Chl2C2
                ad_cff12=ad_cff12+(0.5_r8+adfac)*ad_Chl2C2
                ad_Chl2C2=0.0_r8

!>              tl_cff12=(tl_Bio(i,k,iChlo2)-cff12*tl_fac12)/(fac12+eps)
                adfac=ad_cff12/(fac12+eps)
                ad_fac12=ad_fac12-cff12*adfac
                ad_Bio(i,k,iChlo2)=ad_Bio(i,k,iChlo2)+adfac


!>              tl_fac12=tl_Bio(i,k,iPhyt2)*cff+Bio1(i,k,iPhyt2)*tl_cff
                ad_cff=ad_cff+Bio1(i,k,iPhyt2)*ad_fac12
                ad_Bio(i,k,iPhyt2)=ad_Bio(i,k,iPhyt2)+cff*ad_fac12
                ad_fac12=0.0_r8

!>              tlfac=SIGN(0.5_r8,Chl2C_m(ng)-cff1)
!>              tl_Chl2C=(0.5_r8+tlfac)*tl_cff1+                        &
!>   &                   (0.5_r8-tlfac)*tl_Chl2C_m
                adfac=SIGN(0.5_r8,Chl2C_m(ng)-cff1)
                ad_Chl2C_m=ad_Chl2C_m+(0.5_r8-adfac)*ad_Chl2C1
                ad_cff1=ad_cff1+(0.5_r8+adfac)*ad_Chl2C1
                ad_Chl2C1=0.0_r8

!>              tl_cff1=(tl_Bio(i,k,iChlo1)-cff1*tl_fac1)/(fac1+eps)
                adfac=ad_cff1/(fac1+eps)
                ad_fac1=ad_fac1-cff1*adfac
                ad_Bio(i,k,iChlo1)=ad_Bio(i,k,iChlo1)+adfac

                ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)+0.5*                &
      &                           (ad_cff1/(fac1+eps)+                  &
      &                            ad_cff12/(fac12+eps))
                ad_cff1=0.0_r8
                ad_cff12=0.0_r8

!>              tl_fac1=tl_Bio(i,k,iPhyt)*cff+Bio1(i,k,iPhyt)*tl_cff
                ad_cff=ad_cff+Bio1(i,k,iPhyt)*ad_fac1
                ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+cff*ad_fac1
                ad_fac1=0.0_r8
# else /*2*/

!>              tlfac=SIGN(0.5_r8,Chl2C_m(ng)-cff1)
!>              tl_Chl2C=(0.5_r8+tlfac)*tl_cff1+                        &
!>   &                   (0.5_r8-tlfac)*tl_Chl2C_m
                adfac=SIGN(0.5_r8,Chl2C_m(ng)-cff1)
                ad_Chl2C_m=ad_Chl2C_m+(0.5_r8-adfac)*ad_Chl2C
                ad_cff1=ad_cff1+(0.5_r8+adfac)*ad_Chl2C
                ad_Chl2C=0.0_r8

!>              tl_cff1=(tl_Bio(i,k,iChlo)-cff1*tl_fac1)/(fac1+eps)
                adfac=ad_cff1/(fac1+eps)
                ad_fac1=ad_fac1-cff1*adfac
                ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)+adfac
                ad_cff1=0.0_r8

!>              tl_fac1=tl_Bio(i,k,iPhyt)*cff+Bio1(i,k,iPhyt)*tl_cff
                ad_cff=ad_cff+Bio1(i,k,iPhyt)*ad_fac1
                ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+cff*ad_fac1
                ad_fac1=0.0_r8
# endif /*2*/

!>              tl_cff=tl_PhyCN*12.0_r8
                ad_PhyCN=ad_PhyCN+ad_cff*12.0_r8
                ad_cff=0.0_r8
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!               ========================================================
#ifdef PHYT2 /*1*/
                fac1=AttSW(ng)+AttChl(ng)*Bio1(i,k,iChlo1)+AttFac
                fac12=AttSW(ng)+AttChl(ng)*Bio1(i,k,iChlo2)+AttFac
                fac22=Itop2*(1.0_r8-ExpAtt2)
#else
                fac1=AttSW(ng)+AttChl(ng)*Bio1(i,k,iChlo)+AttFac
#endif
                fac2=Itop*(1.0_r8-ExpAtt)
!               ========================================================
!
#ifdef PHYT2
!>              tl_PAR2=(tl_fac22-PAR2*tl_Att2)/Att2
                adfac=ad_PAR2/Att2
                ad_Att2=ad_Att2-PAR1*adfac
                ad_fac2=ad_fac2+adfac
                ad_PAR=0.0_r8
!
!>              tl_fac2=tl_Itop*(1.0_r8-ExpAtt)-Itop*tl_ExpAtt
                ad_ExpAtt=ad_ExpAtt-Itop*ad_fac2
                ad_Itop=ad_Itop+(1.0_r8-ExpAtt)*ad_fac2
                ad_fac2=0.0_r8
!
!>              tl_Itop=tl_PAR
                ad_PAR=ad_PAR+ad_Itop
                ad_Itop=0.0_r8
!
!>              tl_ExpAtt=-ExpAtt*tl_Att
                ad_Att=ad_Att-ExpAtt*ad_ExpAtt
                ad_ExpAtt=0.0_r8
!
!>              tl_Att=tl_fac1*(z_w(i,j,k)-z_w(i,j,k-1))+               &
!>   &                 fac1*(tl_z_w(i,j,k)-tl_z_w(i,j,k-1))
                adfac1=fac1*ad_Att

!>              tl_PAR=(tl_fac2-PAR*tl_Att)/Att
                adfac=ad_PAR/Att
                ad_Att=ad_Att-PAR1*adfac
                ad_fac2=ad_fac2+adfac
                ad_PAR=0.0_r8
!
!>              tl_fac2=tl_Itop*(1.0_r8-ExpAtt)-Itop*tl_ExpAtt
                ad_ExpAtt=ad_ExpAtt-Itop*ad_fac2
                ad_Itop=ad_Itop+(1.0_r8-ExpAtt)*ad_fac2
                ad_fac2=0.0_r8
!
!>              tl_Itop=tl_PAR
                ad_PAR=ad_PAR+ad_Itop
                ad_Itop=0.0_r8
!
!>              tl_ExpAtt=-ExpAtt*tl_Att
                ad_Att=ad_Att-ExpAtt*ad_ExpAtt
                ad_ExpAtt=0.0_r8
!
!>              tl_Att=tl_fac1*(z_w(i,j,k)-z_w(i,j,k-1))+               &
!>   &                 fac1*(tl_z_w(i,j,k)-tl_z_w(i,j,k-1))
                adfac1=fac1*ad_Att


#else
!>              tl_PAR=(tl_fac2-PAR*tl_Att)/Att
                adfac=ad_PAR/Att
                ad_Att=ad_Att-PAR1*adfac
                ad_fac2=ad_fac2+adfac
                ad_PAR=0.0_r8
!
!>              tl_fac2=tl_Itop*(1.0_r8-ExpAtt)-Itop*tl_ExpAtt
                ad_ExpAtt=ad_ExpAtt-Itop*ad_fac2
                ad_Itop=ad_Itop+(1.0_r8-ExpAtt)*ad_fac2
                ad_fac2=0.0_r8
!
!>              tl_Itop=tl_PAR
                ad_PAR=ad_PAR+ad_Itop
                ad_Itop=0.0_r8
!
!>              tl_ExpAtt=-ExpAtt*tl_Att
                ad_Att=ad_Att-ExpAtt*ad_ExpAtt
                ad_ExpAtt=0.0_r8
!
!>              tl_Att=tl_fac1*(z_w(i,j,k)-z_w(i,j,k-1))+               &
!>   &                 fac1*(tl_z_w(i,j,k)-tl_z_w(i,j,k-1))
                adfac1=fac1*ad_Att
#endif

#ifndef UV_FIXED_TL
                ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-adfac1
                ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+adfac1
#endif
                ad_fac1=ad_fac1+(z_w(i,j,k)-z_w(i,j,k-1))*ad_Att
                ad_Att=0.0_r8

!
!>              tl_fac1=tl_AttSW+tl_AttChl*Bio1(i,k,iChlo)+             &
!>   &                           AttChl(ng)*tl_Bio(i,k,iChlo)
                ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)+AttChl(ng)*ad_fac1
                ad_AttChl=ad_AttChl+Bio1(i,k,iChlo)*ad_fac1
                ad_AttSW=ad_AttSW+ad_fac1
                ad_fac1=0.0_r8
              END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
            ELSE
!             ==========================================================
              fac1=dtdays*NitriR(ng)
              DO k=1,N(ng)
#if defined OXYGEN
                fac2=MAX(Bio2(i,k,iOxyg),0.0_r8)
                fac3=fac2/(K_Nitri(ng)+fac2)
# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                cff3=fac1*fac3*cff
# else
                cff3=fac1*fac3
# endif
#else
                cff3=fac1
#endif
!>              Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
!>              Bio2(i,k,iNO3_)=Bio(i,k,iNO3_)
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
!>              Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
!>              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
!               ========================================================
#ifdef OXYGEN
!>              Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
!>              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                    &
!>   &                            2.0_r8*tl_N_Flux_Nitrifi
                ad_N_Flux_Nitrifi=ad_N_Flux_Nitrifi-                    &
     &                            2.0_r8*ad_Bio(i,k,iOxyg)
#endif
!>              Bio2(i,k,iNO3_)=Bio(i,k,iNO3_)
!>              tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+tl_N_Flux_Nitrifi
                ad_N_Flux_Nitrifi=ad_N_Flux_Nitrifi+ad_Bio(i,k,iNO3_)

!>              tl_N_Flux_Nitrifi=tl_Bio(i,k,iNH4_)*cff3+               &
!>   &                            Bio(i,k,iNH4_)*tl_cff3
                ad_cff3=ad_cff3+Bio(i,k,iNH4_)*ad_N_Flux_Nitrifi
                ad_Bio(i,k,iNH4_)=ad_Bio(i,k,iNH4_)+                    &
     &                            ad_N_Flux_Nitrifi*cff3
                ad_N_Flux_Nitrifi=0.0_r8

!>              Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>              tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio(i,k,iNH4_)*    &
!>   &                             tl_cff3)/(1.0_r8+cff3)
                adfac=ad_Bio(i,k,iNH4_)/(1.0_r8+cff3)
                ad_cff3=ad_cff3+Bio(i,k,iNH4_)*adfac
                ad_Bio(i,k,iNH4_)=adfac
#if defined OXYGEN
# ifdef TDEPENDANCE
!!              tl_cff3=fac1*tl_fac3*cff+                               &
!!   &                  fac1*fac3*tl_cff
!>              tl_cff3=tl_fac1*fac3*cff+                               &
!>   &                  fac1*tl_fac3*cff+                               &
!>   &                  fac1*fac3*tl_cff
                ad_cff=ad_cff+fac1*fac3*ad_cff3
                ad_fac3=ad_fac3+fac1*cff*ad_cff3
                ad_fac1=ad_fac1+fac3*cff*ad_cff3
                ad_cff3=0.0_r8

!>              tl_cff=cff*tl_Bio(i,k,itemp)*LOG(NitriR_t(ng))
#  ifndef UV_FIXED_TL
                ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+                    &
     &                            cff*ad_cff*LOG(NitriR_t(ng))
#  endif
                ad_cff=0.0_r8
# else
!>              tl_cff3=tl_fac1*fac3+fac1*tl_fac3
                ad_fac3=ad_fac3+fac1*ad_cff3
                ad_fac1=ad_fac1+fac3*ad_cff3
                ad_cff3=0.0_r8
# endif
!!              tl_fac3=-fac3/(K_Nitri(ng)+fac2)*tl_fac2
!>              tl_fac3=(tl_fac2-fac3*(tl_K_Nitri+tl_fac2))/            &
!>   &                  (K_Nitri(ng)+fac2)
                adfac=ad_fac3/(K_Nitri(ng)+fac2)
                ad_fac2=ad_fac2-fac3*adfac
                ad_K_Nitri=ad_K_Nitri-fac3*adfac
                ad_fac2=ad_fac2+adfac
                ad_fac3=0.0_r8

!>              tl_fac2=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)))*          &
!>   &                  tl_Bio(i,k,iOxyg)
                adfac=0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg))
                ad_Bio(i,k,iOxyg)=ad_Bio(i,k,iOxyg)+adfac*ad_fac2
                ad_fac2=0.0_r8
#else
!>              tl_cff3=tl_fac1
                ad_fac1=ad_fac1+ad_cff3
                ad_cff3=0.0_r8
#endif
              END DO
!>            tl_fac1=dtdays*tl_NitriR
              ad_NitriR=ad_NitriR+dtdays*ad_fac1
              ad_fac1=0.0_r8
            END IF
!>          tl_PAR=tl_PARsur(i)
            ad_PARsur(i)=ad_PARsur(i)+ad_PAR
            ad_PAR=0.0_r8
          END DO
