
c     input.dat : working folder
c     species_Dl.dat  in folder /IVEA/model
c     Pitzer.para    in folder IVEA/model
c     ML.dat : working folder

c     optional input files
c     parameters.dat ! environmental and solid controlling 1parameters
c     input_back.dat ! control open and closed systems for different species
c     b.fion  ! tune diffusion coefficient for NaCL/NaCO3 solution only
c     acid.dat ! replace oxalic acid by another oganic acid
      
C     35: enzyme concentration aw
c     treat oxalic acid as sulfuric acid for the activity coefficients
c     New compared to v1.0t    
c     More species and more solids
c     oxalic acid, sulfuric acid and phosporic acids
     
c     activity coeffcients of HCO3 and CO3-2
c     with kinetic of CO2 -- HCO3 , iskinetic

c     K1 and K2 for CO2 in pure water is usend sed

c     compared to v2.1b
c     take HCO3 -2 activity coefficient as HSO4- and SO4-2
c     define isupdate ina fcnnew and calhnew
c
c     for the electrodynamic balance (EDB) simulation (imode_output=2) + input_ff.dat is
c     An example is provided in folder RH50_1um_v2.0
c     !!!!!!!!!!!!!
c     flux(I, NP+2): flux from present shell to ammonium oxalate crystal (NP-1)
c     flux(I, lNP-1): flux from inner shellope to ammonium oxalate crystal (NP-1)c     !!!!!!!!!!!!!1
c      
c     input.dat
c line 1:  comment     
C line 2:  0.200000E_b-03 radius in cm
c line 3: 1000		Pressure in hpa  2   nacl f           
c line 4: 5 , 0, 4.35 4 1
c          number of shells, 
c          eddy diffusion coefficient in cm2/s
c     efflorescence Saturation Ratio of solids


c line 5: Imode_output, Idiff, imode_shell 
c         output mode
c         0: (log time), pH , initialize with ML.dat 
c         1: linear output time linear time
c         2: , no pH calculation initialize with ML.dat 
c         3: for background, req
c         5: as 1, but first output time 1E-5 s
c         10+ 0-5, as 0-5, but solid only in the center (ES&T version)
c     22, special treament of NH3NO3 (EDB simulatio), the vapour pressure of NH3 and HNO3 is reduce ny the factor defined in b.para
      
      
c     idiff:
c     0: SLF (modified with b_diff.dat)
c     10: NaCl_sucrose (Aline)
c      13: Na/CO3/Cl- (Bristol)
c     13: NH4NO3/sucrose (Uli)
c     imode_shell
c     0: fixed number of shell
c     1: variable number of shells, call split and merge

      
      
c     line 6: pHmode, iskinetic

c     vmode for pH calculation, 
c     0: default , precise calculation of pH
c     1: many buffers are present
c     2: for the EDB modelling the pH is assumed to be constant (pH=7) 
      
c     iskinetic: 0: fast eqilibrium CO2 + H2O -> H+ HCO3-
c                1: slow kinetic of CO2 + H2O


c     line 7 to 7+N
c          time in seconds
c          relative humidity (RH, 0.5 = 50%)
c          temperature in K
c          acetic acid in the gas phase in ppb
c          NH3 in the gas phase in ppb
c          CO2 in ppm
c          HNO3 in ppb
c          HCl in ppb : > 0 or =0, evaporation considered ; < 0 : no HCl evaporation

c     species.dat
      
      

      implicit real*8 (a-h,m, o-z)
      
      parameter (NSMM=100)
      real* 8 matrix(100,100), matrixb(100),matrixinv(100,100),x13s(100)
      real*8 matrixf(100), xn130(NSMM),dln(nsmm),deltaxn(nsmm)
      real*8   xn13(NSMM) ,a15(NSMM), c15(NSMM)
      real*8 a13(NSMM), c13(NSMM)

      integer iszero (NSMM)
      common /iszero/iszero
      
      
      real*8            xk2tb0(nsmm),xk2t0(nsmm)
      real*8            xk1tb0(nsmm),xk1t0(nsmm)

      real*8 xn1(NSMM),x(NSMM+1),xns(NSMM)
      parameter (np=50,npsolid=12)         ! number of species

            integer lens(Np)      
      integer NSP(npsolid)
      integer NSP_index(npsolid,10)
      real*8 xNSP_nv(npsolid,10)
      COMMON /nsp/ NSP, NSP_INDEX
      COMMON /xnsp/ xNSP_NV

      character*100 path
      common /path/path,Lpath
      
      
      character *60 name(np), filename
      character *200 text
      real*8 cm2(NSMM,NP),velarr(NSMM)
      
      common /cm2/velarr,cm2

      real*8 dl_factor2(2,np)
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin
      
      common  /kout/xfcnnew

      real*8 awshell(NSMM),awshell0(NSMM),
     &     sshell(NSMM,NPsolid),xnsolidshell(nsmm),xjshell(NSMM)
      
    
c      common /awin/ awin
      common /awin/ awin,aws,xvol,xmi
      common/solid/xnsolidt
      common /buffer/  ph3eq,istrace
      
      

      
      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
     & , imode_MA,imode_EDB
      
c       common /isco2/iseqco2,is,nss
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa

      
      common /iskinetic/iskinetic,istrue(NP)
      common /NS0/NS0
      

       common /time/timeA, xn_area
       common /diagnose/ Idiagnose
       common /rcore/ rcore
       
c     species, given by species.dat
      
c     1:  H2O
c     2:  organic 460 molar mass
c     3:  acetic acid total
c     4:  sum of NH3, NH4, and NH4CH3COO
c     5:  sum of CO2 + HCO3-, and CO3--
c     6:  H+
c     7:  OH-
c     8:  CH3COO-
c     9:  undisocciated acetic acid
c     10: NH4CH3COO
c     11: organics which can form protein crystals 
c     12: NH4+
c     13: CO2(aq)
c     14: CO3--
c     15: HCO3-
c     16: Na+
c     17: Cl-
c     18: NO3-
c     19: other cations than H+, Na+, and NH4+ 
c     20: other anions than Cl-, HCO3-, CO3--n , CH3COO-, NO3-
c     21: Virions inactivation with aH+   Influenza virus
c     22: Virions inactivation with aH+   SARS-CoV-2
c     23: NH3(aq)
c     24: H3PO4
c     25: H2PO4-
c     26: HPO4--
c     27: HSO4-
c     28: SO4--
c     29: solid NaCl (the solid shell is put below the liquid, it its assumed that the liquid diffusion is only affected by the presence of the solid due to the size effect ,i.e. the gradient).

      logical ex
       common /sigma/ sigma 


      real*8 flux(NSMM,np+npsolid)



      real*8 ml(NP),ml0(np) ! molality of the species
      real*8 m(NP)  ! molality of the species
      real*4 M4(NP),S4(NPsolid)
      real*8 xn(NSMM,NP),  xnnew(NSMM,NP),  xn0(NSMM,NP)
      real*8 gamma2(NSMM*2,NP),gamma20(NSMM,NP)
      common /gamma2/gamma2
       
c        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
c     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     & , gammah2OA
        real*8 vshellliq(nsmm)

      real*8 ml6shell(NSMM)
      real*8 xhshell(NSMM),phshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell,
     &      vshell0(nsmm)


      real*8 mm(NP) ! molar mass
      real*8 mv(NP) ! molar volume
      integer IZC(NP),ishell(NSMM) !charges of all ions taken from species.dat
      common /m/ mm, mv, izc
      
c     xn1: moles of H2O
c     xn2: moles of sucrose
c     xn3: moles of acetic acid
c     xn4: moles of NH4 
      real*8 ap_eff(NSMM)

      parameter (Ntimes=80000)
      real*8 timetr0(ntimes),rhtr0(ntimes),gash2otr0(ntimes)
      common /ntr0/ Ntr0
      common /ntr/  timetr0,rhtr0
      
      real*8  time1(1)
      real*8 ttr0(ntimes)
      real*8 gasaatr0(ntimes),gasamtr0(ntimes),gasOAtr0(ntimes)
      real*8 gasHNO3tr0(ntimes),gasco2tr0(ntimes),dgasco2tr0(ntimes)
      real*8 gasHCLtr0(ntimes)
      real*8 gaslactr0(ntimes)
      real*8 y1(1),outputtime(ntimes)

      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

      
      real*8 MLL(25000), awL(25000), gcl(25000),gna(25000)
      common /lookup/ MLL,awL,gcl,gna
      real*8 timeEff(100),solidm(100)
      common /enhance/radius,r1,enh_factor
      common /Ienhance/Ienh,iscenter,isAHO

      
      common /timeeq/ timeeq,timeeq0
      common /pi/pi
      common /xn/xn


      real*8  phshell0(nsmm)
      real*8 times
      common /acid/ xkacid1, xkacid2,xhacid,xkacid3

       common /venti/ venti,fheat_11
       common /eutectic/ feutectic
       common /dxmin/ dxmin,factorshell
       save dtime
        common /xm29/ xm29,xlam1_29,xme1_29,xlam2_29,xme2_29
c       common /xvion/ xvion
       common /xvion/ xvion,dl_ion,dl_ions(NSMM)
       common /vap/vap(NP)
       common /ishell/ Ishelleq,imode_shell
       common /fvap/ alpha1,alpha0,ppvap,fvap_factor,vapratio        
       fvap_factor=1d0
       alpha1=1
       
       xvion=195
       ilate=0d0
       
c     averaging for output
       svap4=0d0
       svap18=0d0
       sml12=0d0
       sml18=0d0
       sml2=0d0
       sstime=0d0

       fheat_11=1d0

       
       istimeeq=0
       isdiv=0
       feutectic=-1
       venti=1d0
       xdissolid=5D-4
       fbutter=0.01d0
       npliq=np-npsolid
       npl=np-npsolid
       sdtime=0d0
       dtimesplit=0d0
       dxmin=1D-7
       isolidmin=0
       ph3eq= 4d0
       imin=0
       jmin=0
       DO I=1,nsmm
          do J=1,npsolid
             sshell(I,j)=0d0
             phshell(i)=10d0
          enddo 
          enddo
       
       DO J=1,np
          iszero(j)=0
       enddo
       
c       path='/cluster/home/bluo/'
c       path='/home/bluo/'
       path='/home/medina/'
       do i = 100,1,-1
          if (path(i:i) .ne. ' ') then
             lpath= I
             goto 478
             endif
          enddo
 478      continue
          xn40=-10
          xn180=-10
          
          
          INQUIRE (FILE='droplet.dat', exist=ex)
c          timenuc_drop=-100
          timeeff(1)=-1000d0
          nnuc=0
          
          if (ex) then
             open(99,FILE='droplet.dat')
            read(99,*) fheat_11



c          print*, 'latent heat', fheat_11



          
          
          DO I=1,100
          read(99,*,end=126) ,timeeff(I),solidm(I)
        enddo

 126    nnuc=I-1
        
        Inuc=1
        
          
          close(99)
       endif
          
        Inuc=1


       istrace=1                ! reduce the flux by istrace
c     1: no reduction
       
       INQUIRE (FILE='parameters.dat', exist=ex)
       
       if (ex)then
       open(99,file='parameters.dat')
       read(99, *) feutectic    ! thickness of eutectic structure in cm
c     ! <  0d0 , no eutektikum corretion
       if ( feutectic.le. 1D-7 .and. feutectic .gt.0d0) then
          print*,
     &'error:THickness of eutectic lamells less than 1nm, stop '
          print*, feutectic
          stop
       endif
       
       read(99,*) venti         ! gas phase ventilation factor
c     ! <  0d0 , no ventivation correction
       if (venti.le.1d0) venti=1d0
c     distance between solid shells, 
         read(99,*) xdissolid
      if ( xdissolid.le.0d0) xdissolid=5D-4
      if ( xdissolid.gt.0d0 .and. xdissolid.le.10E-7) then

         print*, 'Distance < of solid shells <  50nm'
         print*, 'error: input a value > 50E-7 cm stop'

         print*, xdissolid
         stop
         
      endif
      factorshell=10d0
        read(99,'(A)')text
        read(text,*)         dxmin
        read(text,*,end=231,err=231)  dxmin, factorshell
 231    continue
        print*, dxmin, factorshell

        

      
         Idiagnose=0
         read(99,*,end=991) Idiagnose

         read(99,*,end=991) ph3eq
c         read(99,*,end=991) istrace
         read(99,*,end=991) xvion
         read(99,*,end=991) alpha1
         alpha0=1d0
         read(99,*,end=991) alpha0
         ppvap0=0d0
         
         read(99,*,end=991) ppvap0 ! pNH4.pHNO3 ppb**2

         ppvap1=0d0
         
         read(99,*,end=991) ppvap1 ! pNH4.pHNO3 ppb**2
         dppvap=0
         read(99,*,end=991) dppvap ! pNH4.pHNO3 ppb**2
         
         taup=1d0
         read(99,*,end=991) taup ! pNH4.pHNO3 ppb**2
         
         vapratio=1
         read(99,*,end=991) vapratio 
         fvap_factor=1d0
         read(99,*,end=991) fvap_factor

         
 991     close(99)
       endif
       
       

       

 429  continue

      dxmin0=dxmin
        isolidmin=0      
      iseq3 =0
      iseqnh3=0
      N1=1
      iscenter=0
      solidsm=0
      Nxhh=1
      loop=0

      loopkka=0
      Ieffcycle=1
      Is=0
      DO I=1,2*nsmm
         DO J=1,np
      gamma2(i,j)=1d0
      enddo
      enddo
      
      xloop=0
      
      
c     IS the shell number for calhnew and calhnewCO2
c     if IS .eq. 0  or MH+ <= 0 set the the default range for H+ in dbrent  
c     if IS  > 0, then set the range to MH+/10 to MH+*10 to save CPU time
c     this is still not full tested


      sigma = 72                ! surface tension erg/cm2 of water
      aweff = 0.569  ! efflorescence water activity (not valid for output mode 2 and 12), for one can moddify this value, it NaCl nucleate at different RH
      pi=dacos(-1d0)
      xnsolid=0d0
c      xdissolid=5E-4            !  we set the miniumum distance between shells with solid NaCl be 5 um here
c     If the NaCl dentrite structure is larger or smaller than 5 um given here, one can modify this value here.


c     ------------------------------------------------------------------------------------------------------------------------------------------------      
c     Read input data
c     ------------------------------------------------------------------------------------------------------------------------------------------------      

      DO i= 1, 130
         phh=.1d0*i
               xmm=0
               dd= tau_ivea(phh, xmm)

               write(62,*) phh,dd*dlog(100d0)

               call cal_tau_sars(T,RH,phH,taus)

       write(63,*) phh,taus*dlog(100d0)
       
             enddo

             filename='species_50.dat'
             open(1,file=filename)
             
             read(1,'(A)') text 

      DO I=1,NP
             read(1,'(A)') text 
             do j=1,30
             if (text(J:J) .eq. ',') goto 591
             enddo
          print*, 
     & 'error: after species name, a comma should follow!'

          stop

 591         filename=text(1:J-1) 
              text= text(j+1:200)

          read(text,*) ii , MM(i),
     & mv(i), dl_factor2(1,i) ,izc(i)
          dl_factor2(2,i)=dl_factor2(1,i) ! 
          if (I.gt. np-npsolid) then
             ip = I-np+npsolid
c             print*,ip
             
             read(text,*) ii , MM(i),
     & mv(i), dl_factor2(1,i) ,izc(i), nsp(Ip) ! 
c             print*,ip,nsp(ip)
             
              
              read(text,*) ii , MM(i),
     &             mv(i), dl_factor2(1,i) ,izc(i), iii
     & , (nsp_index(Ip,k),xnsp_nv(Ip,k), k=1,nsp(ip))
          endif
          if (I.le.npliq) then
          write(6,'(A20,I5, 3F15.4,I5)') filename,i,MM(i),
     & mv(i), dl_factor2(1,i) ,izc(i)
       else
       write(6,'(A20,I5, 3F15.4,2I5,10(I5,F5.1) )') filename,i,MM(i),
     &         mv(i), dl_factor2(1,i) ,izc(i), nsp(Ip)
     & , (nsp_index(Ip,k),xnsp_nv(Ip,k), k=1,nsp(ip))
          

       endif

          name(I) = (filename)
c     remove '.' in name
          DO j=1,40
             if(filename(j:j).eq.'.') then
             
             filename(j:40)=filename(j+1:40)
          endif
          enddo
          lens(I)=len((name(I)))
          
          print*, 'A'//filename(1:lens(I))//'B', lens(I)
          
          name(I)=filename
          
          
             enddo
             close(1)

c     relative to H2O
             ff=dl_factor2(1,1)
             DO I=1, Np
         dl_factor2(1,i)= dl_factor2(1,i)/ ff
c         print*, i, dl_factor2(1,i)
         write(27,*) I, izc(I)
      enddo
                

             
c     read the parameter for the replacement of the oxalic acid
c     MM, MV, H*, xk1, xk2
      INQUIRE (FILE='acid.dat', exist=ex)

      xkacid1=-1
      xkacid2=100
      xkacid3=100
      xhacid=-1

         xlam1_29=0d0
         xme1_29=1d0
         xlam2_29=0d0
         xme2_29=5D0
      if (ex) then
         open(1,file='acid.dat')
        read(1,'(A)')text
        read(text,*)         xmm,xmv,xkacid1, xkacid2,xhacid
         xkacid3=100

        read(text,*,end=823)xmm,xmv,xkacid1, xkacid2,xhacid,xkacid3
 823    xkacid1=10.d0**(-xkacid1)
         xkacid2=10.d0**(-xkacid2)
         xkacid3=10.d0**(-xkacid3)

        if(xkacid2.le.1D-50) xkacid2=0d0
        if(xkacid3.le.1D-50) xkacid3=0d0
    
         
         

         read(1,*,end=124) ii,xlam1_29
         read(1,*,end=124) ii,xme1_29
         read(1,*,end=124) ii,xlam2_29
         read(1,*,end=124) ii,xme2_29
         
         
         
 124      continue
         print*,'xme2_29', xme2_29
         print*,'xme1_29', xme1_29
         

          
         
         mm(29)= xmm
         mm(30)= xmm-mm(6)
         mm(31)= xmm-2*mm(6)

         mv(29)= xmv
         mv(30)= xmv-mv(6)
         mv(31)= xmv-2*mv(6)
         print*, mm(29)
         

         
 123     continue
         close(1)
      endif


      

             open(1,file='input.dat')
      read(1,'(A)') filename

             if (filename(1:2).eq.'pH') then
                iskinetic=0
                IS=1
           read(filename(4:20),*)        T
           print*,t
           close(1)
           open(1,file='ML.dat')
           ml=0d0
           DO J=1,np
              read(1,*) ii,ml(j)
           enddo
              call calhnew(t,ml)
              DO J=1,np
          if (ml(j).gt.0d0) write(6,'(A20,I5,10E15.6)')name(j), j, ml(j)
                 enddo
      call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      pH = dlog( ml(6)*gammaH)/dlog(.1d0)
      print*, ' pH ', pH
      print*, ' aw ', aw

      
      INQUIRE (FILE='buffer.dat', exist=ex)

      
      if (ex) then 
      open(2,file='output_pH.dat')
      open(10,file='buffer.dat')
      
      read(10,'(A10)')
      
      close(3)
      close(3)
      dO ii=1,100,1
         read(10,*,end=389) ph0, xnacl, xca1,xca2, xnap1,xnap2
         ml=0d0
         xvol=  xca1+xca2+ xnap1+xnap2

         xmp = xnap1*0.2/1000+ xnap2*0.2/1000 ! mols P
         xvp = xnap1*0.2/1000*mv(26) + xnap2*0.2/1000 *mv(37) !cm3
         xmNa =  2*xnap1*0.2/1000+ 3*xnap2*0.2/1000 ! mols Na
         xvNa = xmNa*mv(16) 
         xvpw = xnap1 + xnap2 - xvp-xvNa

         xmc = (0.1*xca1+xca2)/1000 ! mols citric acid
         xvc = (0.1*xca1+xca2)/1000*mv(29) ! cm3 volume citric acid
         xvcw = xca1+xca2 - xvc

         xvw  =xvcw+xvpw

         xmw=xvw * mm(1)/mv(1)
         xmnacl=0d0
         if (xnacl.ge.2)  xmnacl = 3.945 /(mm(17)+mm(16))

         xmnacl = xmnacl * 1000/xmw
         ml(24)= xmp*1000/xmw
         ml(29)= xmc*1000/xmw
         ml(16)=xmNa*1000/xmw
         
         print*, 'vol', xvol, xvw
         
         ml(17) =xmnacl
         ml(16) =ml(16)+xmnacl
         
         ml(1)=1000/mm(1)
         
         
c     xvol=1000d0
c         xvol=1
         
c         ml(29)= (0.1*xca1+xca2) /xvol
c         ml(24) =(0.2* xnap1+0.2*xnap2)/xvol
c         xna=(0.4* xnap1+0.6*xnap2)/xvol
c         xnacl = xnacl*(1000-) /1000d0

c         ml(17)= xnacl
c         ml(16)= xnacl + (0.4* xnap1+0.6*xnap2)/xvol
c         vv= 1000- mv(16)*ml(16)- mv(17)*ml(17)- mv(24)*ml(24)
c     & - mv(29)*ml(29)         


         

c         ml = ml* 1000/vv       ! convert Molarity to molality
c        ml(1)=1000/mm(1)
         ml0=ml
         
      call calhnew(t,ml)
      call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      pH = dlog( ml(6)*gammaH)/dlog(.1d0)


      
      write(6,'(I6,4F15.6,10E15.6)')Ii,ph,ml0(16),ml0(17),
     &ml0(24),ml0(29),
     &     ml(24)+ml(25)+ml(26)+ml(37),ml(24),ml(25),ml(26),ml(37)
     & ,gamma2(1,6)
      write(2,'(I6,4F15.6,10E15.6)')Ii,ph,ml0(16),ml0(17),
     &ml0(24),ml0(29),
     & ml(24)+ml(25)+ml(26)+ml(37),ml(24),ml(25),ml(26),ml(37)
     & ,gamma2(1,6)

      write(3,'(I6,4F15.6,10E15.6)')Ii,ml0(16),ml0(17),
     &ml0(24),ml0(29), ph0

      
      enddo
           
 389   print*, ' Output is in output_pH.dat'
       call cpu_time(time11)
       DO xx=1,1D3,1
      call calhnew(t,ml)
      call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
       enddo
       call cpu_time(time2)
       print*,time2-time11
      endif
      
                stop
      
                
                stop
             endif


      read(1,*) r0
      read(1,*) press0
      read(1,*) NSmax,Deddy, (AP_eff(JJ),jj=1,npsolid)
      ns0=nsmax
      

c     for (NH4)2C2O4


      
      isNH3close=0
      isHNO3close=0
      isOAclose=0
      isaaclose=0
      isHCLclose=0
      isLACclose=0
      
      pm=0

      
      INQUIRE (FILE='input_back.dat', exist=ex)
      if (ex) then
      open(20,file='input_back.dat')
         read(20,*,end=1)  pm ! particular matter in cm3/m3
         read(20,*)  isNH3close ! bei 50 ppbHNO3, is NH3 close=1
         read(20,*)  isHNO3close ! bei 50 ppbHNO3, is NH3 close=1
         read(20,*)  isOAclose ! bei 50 ppbHNO3, is NH3 close=1, final size
         read(20,*)  isHClclose ! bei 50 ppbHNO3, is NH3 close=1, final size
         read(20,*)  isAAclose ! bei 50 ppbHNO3, is NH3 close=1, final size
         read(20,*)  islacclose ! bei 50 ppbHNO3, is NH3 close=1, final size
 1       print*, 'is close 1', isnh3close, isHNO3close
         if(pm.le.0d0) then
      isNH3close=0
      isHNO3close=0
      isOAclose=0
      isaaclose=0
      isHCLclose=0
      endif
         endif

      
        isclose= isNH3close+ isHNO3close+ isOAclose+ isAAclose+
     & isHCLclose+ islacclose

        idiff=0
        imode_shell =0
        imode_NH4NO3 =0
        
        read(1,'(A)')text
        read(text,*)         Imode_output, idiff
        read(text,*,end=131)         Imode_output, idiff   , imode_shell           ! output mode , 0: log(time), 1: linear time
 131    continue
        iscenter=0
        if ( imode_output.ge.40) then
           imode_output=imode_output-40
           imode_EDB=1
        endif

        
        if ( imode_output.ge.30) then
           imode_output=imode_output-30
           imode_MA=1
        endif

        
         
        if ( imode_output.ge.20) then
           imode_output=imode_output-20
           imode_NH4NO3=1
        endif

      if ( imode_output.eq.3) aweff=0d0




      if (imode_output.ge.10) then
         iscenter=1
         Imode_output=Imode_output-10
      endif
      
      imode_outputt=      imode_output  
      if (imode_output.ge.3)  imode_outputt=0
      
      ismode=0
      if (imode_output.eq.0) ismode=1
      if (imode_output.eq.1) ismode=1
      if (imode_output.eq.2) ismode=1
      if (imode_output.eq.3) ismode=1
      if (imode_output.eq.4) ismode=1
      if (imode_output.eq.5) ismode=1
      if (ismode.eq.0) then
         print*,
     & 'error: Input a valid output mode (0,1,2,4,5,10,11,12,14,15)'
         
         stop
      endif
      


      
 26     continue

c     fenzyme
        fCl=0d0
c     catkm:  lg_10(kcal/km) ph 7
c     dcatkm  dlg(kcal/km)/dpH   1/pH
        
        read(1,*) Imode_ph,iskinetic, catkm, dcatkm, fCO3
        print*, 'iskinetic', iskinetic
        imode_eq=0
        if (imode_ph.ge.10) then
           imode_eq=1
           imode_ph=imode_ph-10
        endif
        
      if (imode_ph.lt.0 .or. imode_pH.gt.3) then
         print*, 'Input a valid pH mode (0,1,2) '
         stop
      endif

c     ------------------------------------------------------------------------------------------------------------------------------------------------   c     Idiff ist set here to 0   
c     For sensitivity calculations, one has the change the Idiff here
c     ------------------------------------------------------------------------------------------------------------------------------------------------


c      
c line 6: 0
c         diffusion mode Idiff
c         0: activity gradient, with charge balance except (H+,OH-), SLF diffusivity
c         1: activity gradient,       SLF diffusivity, w.o. charge balance
c         2: concentration gradient,  SLF diffusivity, w.o. charge balance
c         3: concentration gradient,  sucrose diffusivity, with charge balance
c         4: concentration gradient,  citric acid diffusivity, w.o. charge balance.
c         5: Dl_ions = D_H2O_walker     activity with charge balance 
c         6: Dl_ions = D_H2O_thiswork   activity with charge balance


        Iexit=1 ! the index when the time=0, the exhaled air exits the nose or mouth
c     For the mixing, the air mixed air at time=0 with the ambient air.
        

        DO I=1,ntimes
         read(1,*,end=222) timetr0(I), rhtr0(I),ttr0(I),
     &  gasaatr0(I),gasamtr0(I),gasco2tr0(I),gashno3tr0(I),gashCLtr0(I)
     & ,gasOAtr0(I) ,gaslactr0(I)
         T=ttr0(I)
         gash2otr0(I) = vwater(t)*rhtr0(I)
         if (dabs(timetr0(I)).le.1D-10) Iexit=I
c         print*,i, timetr0(I), ttr0(I)

      enddo
 222     Ntr0=I-1
         print*, ntr0
c time derivative of gas CO2         

         
         dgasco2tr0(1)=0d0
         DO I=2,Ntr0
            dgasco2tr0(I)= gasco2tr0(I)- gasco2tr0(I-1)
         enddo

         close(1)
       
       
c     test diffusion coefficients in the liquid for EDB simulation

cccccccccccccccccccccccccccccccccccccaccccccccc      
         
       
         
         if (imode_output.eq.2) then

            t=293.15
        DO i =0,100, 1
           aw= i*.01
           call cal_dlaw(T,aw,dl)
       ddcl=dl_factor2(1,17)
       dd=dl_factor2(1,16)
       call cal_dlaw_walker(T,aw,dlw)
       call cal_dlaw_suc(T,aw,dlsuc)
       call cal_dlaw_citric(T,aw,dlcitric)
       call cal_dlaw_walker_h2o(T,aw,dlm)
c       write(62,'(15E15.6)')aw,dlm, dl*dd, dl*ddcl,dlcitric,DLSuc,dlw
      enddo
      endif

      open(34,file='output_n.dat')
      is=1
               DO I=1,100
               ML=0
            ml(1)=1000d0/mm(1)
            xmm=1D-5+ I
                  ml(2)=0d0
                  if( imode_NH4NO3.eq.1) then
               ml(4)= xmm
               ml(12)= xmm
               ml(18)= xmm
            else

               ml(29)=xmm

            endif
            look=1
            tt=298
              call calhnew(tt,ml)
      call aw_back
     & (tt,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      write(6,'(A,6E15.6)'), 'L',xmm, ml(30),aw, gammah,gammaoa,gammahoa
              call calhnew_model(tt,ml)
      call aw_back_model
     & (tt,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      write(6,'(A,6E15.6)'), 'M',xmm, ml(30),aw, gammah,gammaoa,gammahoa

      

      
      enddo
      
c      stop
      
      if (idiff.eq.13) then


         open(2,file='ML.dat')
         DO I=1,NP
            read(2,*) II,ml(I)
            if (ML(I).le.1D-40) ml(I)=0d0
         enddo
         close(2)
      endif

      open(21,file='output_shells_Dl.dat')
      
      open(75,file='output_Dl.dat')
       print*,'befor dl'
       DO ii=0, 100
          aw=0.01*ii

         call cal_h2o_guess(T,aw,dlH2O)
         write(24,*) aw, dlH2O


       xvv=xvion
       call cal_dlaw_eta(T,aw,dl,xvv)
          
c       call cal_dlaw(T,aw,dl)  !diffusin coefficient of ions
c       write(64,'(15E15.6)') aw,dl
       xvv=mv(1)
       call cal_dlaw_eta(T,aw,dlm,xvv)

c          call cal_dlaw_walker_h2o(T,aw,dlm)
       xvv=mv(29)
c       xvv=368
       call cal_dlaw_eta(T,aw,dlacid,xvv)

       xvv=mv(2)
       call cal_dlaw_eta(T,aw,dlmatrix,xvv)
       
c       write(65,'(15E15.6)') aw,dlm,dl,dlacid,dlmatrix

       if (idiff.eq.13) then

          xmv= mv(1)
        call  cal_dl_an_suc(T,ML,aw,xmv,dlm)    
          xmv= xvion
          call cal_dl_an_suc(T,ml,aw,xmv,dl)      
          xmv= mv(2)
          call cal_dl_an_suc(T,ml,aw,xmv,dlmatrix)      

          xmv= 368              ! mv(29)
          call cal_dl_an_suc(T,ml,aw,xmv,dlacid)      
          
       endif

        call cal_dlaw(T,aw,dl0)  !diffusin coefficient of ions
c       call cal_dlaw_walker_mod(T,aw,dlm0)
        call cal_dlaw_walker_mod(T,aw,dlm0)  ! present work
       
       
       write(75,'(15E15.6)') aw,dlm,dl,dlacid,dlmatrix,t,Dlm0,dl0
       
       enddo
       print*,'finish dl'
       DO ww =0d0,1d0,0.01d0
          write(50,*) ww, alpha0+ (alpha1-alpha0)*ww
       enddo

       
c     sucrose
       if (Idiff.eq.13) then
       open(75,file='output_Dl_suc.dat')
       DO ii=0, 100
          aw=0.01*ii
          Ml=0
          ML(1)=1000/MM(1)
          ML(2)=1D0


          xmv= mv(1)
        call  cal_dl_an_suc(T,ML,aw,xmv,dlm)    
          xmv= xvion
          call cal_dl_an_suc(T,ml,aw,xmv,dl)      
          xmv= mv(2)
          call cal_dl_an_suc(T,ml,aw,xmv,dlmatrix)      
          xmv= mv(29)
          call cal_dl_an_suc(T,ml,aw,xmv,dlacid)      
          

       
       write(75,'(15E15.6)') aw,dlm,dl,dlacid,dlmatrix
       enddo

      endif
      
       if (idiff.eq.13) then
          IS=1
          open(75,file='output_vap_MA.dat')
          open(76,file='output_mfs.dat')
          open(77,file='output_mfs_0.7.dat')
       T1=298.15
       
c     pure NH4NO3
c       T=298.16
       DO xxx=.01, 100,.1
          ml=0
          ml(1)=1000/mm(1)
          ml(4) = xxx
          ml(12) = xxx
          ml(18) = xxx
c          xmfs= ml(29)*mm(29)/(1000+ ml(29)*mm(29))
       call vapnew(T1,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa1,pl)
c       write(16,'(7E15.6)') aw, pnh3*pHNO3,xxx,xxx
       
       
       

       
          ml=0
          ml(1)=1000/mm(1)
          ml(29) = xxx
          xmfs= ml(29)*mm(29)/(1000+ ml(29)*mm(29))
       call vapnew(T1,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa1,pl)
       write(76,'(7E15.6)') aw, xmfs,gamma2(1,29)

       
c       call vapnew(T1,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa1,pl)
          ml=0
          ml(1)=1000/mm(1)
          ml(29) = xxx*0.7d0
          ml(2) = xxx*0.3d0
          xms= ml(29)*mm(29)+ml(2)*mm(2)
          xmfs= xms/(1000+ xms)

       write(77,'(7E15.6)') aw, xmfs,gamma2(1,29)
       

       
       ml=0
       ratio=0.5
          ml(1)=1000/mm(1)
          ml(29) = xxx*ratio
          ml(2)= xxx*(1-raio)

          call vapnew(T1,ML,aw5,pacetic,pnh3,pHNO3,PHCL,PCO2,poa5,pl)
       ml=0
       ratio=0.75
          ml(1)=1000/mm(1)
          ml(29) = xxx*ratio
          ml(2)= xxx*(1-raio)

          call vapnew(T1,ML,aw75,pacetic,pnh3,pHNO3,PHCL,PCO2,poa75,pl)
       ml=0
       ratio=0.25
          ml(1)=1000/mm(1)
          ml(29) = xxx*ratio
          ml(2)= xxx*(1-raio)

          call vapnew(T1,ML,aw25,pacetic,pnh3,pHNO3,PHCL,PCO2,poa25,pl)

          
       
          write(75,'(7E15.6)') aw, poa1, 1d0
          write(75,'(7E15.6)') aw75, poa75,.75d0
          write(75,'(7E15.6)') aw5, poa5,0.5d0
          write(75,'(7E15.6)') aw25, poa25,0.25d0
           
      enddo
      endif
      print*,'aa'
c      stop
      

       if (idiff.eq.13) then

          open(75,file='output_Dl_AN.dat')
       
c     pure NH4NO3
       ML(2)=0d0
       ML(12)=1
       ML(4)=1
       ML(18)=1
       
       DO ii=0, 100
          aw=0.01*ii

       
          xmv= mv(1)
        call  cal_dl_an_suc(T,ML,aw,xmv,dlm)    
          xmv= xvion
          call cal_dl_an_suc(T,ml,aw,xmv,dl)      
          xmv= mv(2)
          call cal_dl_an_suc(T,ml,aw,xmv,dlmatrix)      
          xmv= mv(29)
          call cal_dl_an_suc(T,ml,aw,xmv,dlacid)      
       write(75,'(15E15.6)') aw,dlm,dl,dlacid,dlmatrix
       enddo

       
c     pure sucrose
       ml=0
       ml(1)=1000/mm(1)
       
       open(75,file='output_Dl_H2O.dat')

          DO ii=0, 100
          aw=0.01*ii

          ML(2)=1d0
          ml(4)=0d0
          ml(12)=0d0
          ml(18)=0d0
       
          xmv= mv(1)
        call  cal_dl_an_suc(T,ML,aw,xmv,dl00)    
          xx=.1
          ML(2)=1-xx
          ml(4)=xx
          ml(12)=xx
          ml(18)=xx
        call  cal_dl_an_suc(T,ML,aw,xmv,dl01)    
          xx=.3
          ML(2)=1-xx
          ml(4)=xx
          ml(12)=xx
          ml(18)=xx
        call  cal_dl_an_suc(T,ML,aw,xmv,dl03)    
          xx=.5
          ML(2)=1-xx
          ml(4)=xx
          ml(12)=xx
          ml(18)=xx
        call  cal_dl_an_suc(T,ML,aw,xmv,dl05)    
          xx=.7
          ML(2)=1-xx
          ml(4)=xx
          ml(12)=xx
          ml(18)=xx
        call  cal_dl_an_suc(T,ML,aw,xmv,dl07)    
          xx=1d0
          ML(2)=1-xx
          ml(4)=xx
          ml(12)=xx
          ml(18)=xx
        call  cal_dl_an_suc(T,ML,aw,xmv,dl10)    
       write(75,'(15E15.6)') aw,dl00,dl01,
     & dl03,dl05,dl07,dl10
      enddo

          open(75,file='output_Dl_ions.dat')
          open(76,file='output_Dl_29.dat')

          DO ii=0, 100
          aw=0.01*ii
          ml=0
          ml(1)=1000/mm(1)
          
          ML(2)=1d0
          ml(4)=0d0
          ml(12)=0d0
          ml(18)=0d0
       
          xmv= xvion
        call  cal_dl_an_suc(T,ML,aw,xmv,dl00)    
          ml=0
          ml(1)=1000/mm(1)
          ML(2)=1d0

        xmv= mv(29)
        call  cal_dl_an_suc(T,ML,aw,xmv,dl00_29)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.3
          ML(2)=1-xx

          ml(4)=xx/2
          ml(12)=xx/2
          ml(18)=xx/2
          
          xmv= xvion
        call  cal_dl_an_suc(T,ML,aw,xmv,dl01)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.3
          ML(2)=1-xx
          ml(29)=xx
        xmv= mv(29)

        call  cal_dl_an_suc(T,ML,aw,xmv,dl01_29)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.5
          ML(2)=1-xx
          ml(4)=xx/2
          ml(12)=xx/2
          ml(18)=xx/2
          xmv= xvion
          call  cal_dl_an_suc(T,ML,aw,xmv,dl03)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.5
          ML(2)=1-xx
          ml(29)=xx
        xmv= mv(29)
        call  cal_dl_an_suc(T,ML,aw,xmv,dl03_29)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.7
          ML(2)=1-xx
          ml(4)=xx/2
          ml(12)=xx/2
          ml(18)=xx/2
          xmv= xvion
          call  cal_dl_an_suc(T,ML,aw,xmv,dl05)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.7
          ML(2)=1-xx
          ml(29)=xx
        xmv= mv(29)
        call  cal_dl_an_suc(T,ML,aw,xmv,dl05_29)    

          ml=0
          ml(1)=1000/mm(1)
        
        xx=.8
          ML(2)=1-xx
          ml(4)=xx/2
          ml(12)=xx/2
          ml(18)=xx/2
          xmv= xvion
        call  cal_dl_an_suc(T,ML,aw,xmv,dl07)    
          ml=0
          ml(1)=1000/mm(1)
          xx=.8
          ML(2)=1-xx
          ml(29)=xx
        xmv= mv(29)
        call  cal_dl_an_suc(T,ML,aw,xmv,dl07_29)    

          ml=0
          ml(1)=1000/mm(1)

        xx=1d0
          ML(2)=1-xx
          ml(4)=xx
          ml(12)=xx
          ml(18)=xx
          xmv= xvion
        call  cal_dl_an_suc(T,ML,aw,xmv,dl10)    
          ml=0
          ml(1)=1000/mm(1)
          xx=1
          ML(2)=1-xx
          ml(29)=xx
        xmv= mv(29)
        call  cal_dl_an_suc(T,ML,aw,xmv,dl10_29)    


        write(75,'(15E15.6)') aw,dl00,dl01,
     & dl03,dl05,dl07,dl10,t,xvion
        write(76,'(15E15.6)') aw,dl00_29,dl01_29,
     & dl03_29,dl05_29,dl07_29,dl10_29
      enddo
      
      
      endif
          

       call flush(64)
       call flush(65)
       call flush(75)

cccccccccccccccccccccccccccccccccccccccccccccc      
       if (idiff.eq.13) then
          open(75,file='output_Dl_AN.dat')
      ML(2)=0d0
      
      
       endif
       
      

c     if time > timeeq, it is assumed that the water activity of the droplet is equal to RH
c     no kinetic calculation for H2O

       
        timeeq=1d10
        timeeq3=1D10

c     H2O Equilibrium for exhaled aerosol

        if (imode_output.ne.2) then 

        if (r0.le. .3D-4) then
           timeeq=2
           endif
        if (r0.le. .022D-4) then
           timeeq=1
           endif
           isc=isHNO3close+isHCLclose+isNH3close+isOAclose+isAAclose

           if (r0.ge. 0.3E-4 .or. isc .ge. 1 )then
              timeeq=2* (r0*1D4)**2
              endif
             endif
              
c              timeeq=1D20
              
              timeeq0=timeeq

              if (iskinetic.eq.1.and.NS.eq.1)  timeeq=1D10

              



c     define output times
       time11=1D-2
       if (imode_output.eq.5)        time11=1D-5
       if (imode_output.eq.5)      imode_output=0

       time22= timetr0(NTR0)
        Noutput=1001
       DtimeL = dlog(time22/time11)/(Noutput-1)
      
       


       
      if ( imode_output .eq.0 .or. imode_output.eq.4
     & .or. imode_output.eq.3 )then

         if (TIMETR0(1).ge.0d0) then

         Noutput=1001
         do i=1,nOUTPUT
         outputtime(I)= time11* dexp((I-1)* dtimeL)
         enddo

         else
         Noutput=1101
         dtime= (TIMETR0(2)-TIMETR0(1))/99d0
         

         do i=1,100
          outputtime(I)= timetr0(1)+ dtime*(I-1)
c          print*, i, outputtime(I)
          enddo


         do i=1,nOUTPUT
         outputtime(I+100)= time11* dexp((I-1)* dtimeL)
         enddo
         endif

c         if (i.le.100)print*, I, outputtime(I)

         endif



      if ( imode_output.eq.1   ) then
         noutput= 5001
         
       Dtime = (time22-time11)/(Noutput-1)
         do i =1,noutput
            outputtime(I)= TIMETR0(1) + (i-1)*DTIME
         ENDDO
         
         endif


         
            if (imode_output .eq. 2  ) then
c               timetr0(1)=timetr0(2) -1d0
               Noutput = Ntr0
               DO I=1,Ntr0
             outputtime(I)= timetr0(I)
               enddo
                           endif

             outputtime(Noutput+1)=  outputtime(Noutput)+100
                print*, noutput
          





c     ----------------------------------------------------------------------
c     Initialization  
c     ----------------------------------------------------------------------

                DO I=1,np
                   print*,i,mv(i)
                   enddo

        t0=ttr0(1)
        Ta=T0
        tdrop=T0
        T= t0
        rh=rhtr0(1)
        Amisch = apNaCL(t)
        press=press0

c     read the inital compostion (z.B., given SLF recipe or equilibrated with given gas phases)

        

         open(1,file='ML.dat')
         DO I=1,NP
         read(1,*) ii,ml(I)
c         print*,i,ml(i)
         
      enddo
      t=ttr0(1)
      HCO2 = 0.034 * dexp(2300 *(1/T-1/298.15))
                xm13=hco2/1013.5*press0

      if (gasco2tr0(1).ge.1D-30 ) then
         if (Ml(5).le. 1D-10) ml(5)=1D-10
         
      endif
      

      iskinetic0=iskinetic
      iskinetic=0
      call calhnew(t,ml)
      iskinetic=iskinetic0
      
      
      
          if (Imode_NH4NO3.eq.1) then
           call aw_back
     & (t,ML,aw0,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

         call cal_ml(T,rh,ML)
           call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           print*,aw0,rh,aw
c           print*, rhtr0(1), rhtr0(2)
           
        endif

      
c     if (ML(11).le.1D-40)ML(11)=1D-40 ! 11 take as a trace of inert virus 
c     define which species is zero
c         read(1,*,end=222) timetr0(I), rhtr0(I),ttr0(I),
c     &  gasaatr0(I),gasamtr0(I),gasco2tr0(I),gashno3tr0(I),gashCLtr0(I)
c     & ,gasOAtr0(I) ,gaslactr0(I)
            
c     specical treament of CO2

c     define istrue:1: species > 0
c     define istrue:0: species =< 0
      iskinetic0=iskinetic

      DO J=1,NP
      istrue(J)=1
      if (ml(j).le.1D-40 ) istrue(j)=0
      if (ml(j).le.1D-40 ) ml(j)=0d0
      enddo
      
c      if (ML(2.le.0d0) istrue(2)=0d0
      
c     consider volatile species
      
      gasaa=0d0
      gasAm=0d0
      gasCO2=0d0
      gasHNO3=0d0
      gashcl=0d0
      gasOA=0d0
      gasla=0
      DO I=1, ntr0
         if ( gasaatr0(I).gt. gasaa) gasaa= gasaatr0(I)
         if ( gasamtr0(I).gt. gasam) gasam= gasamtr0(I)
         if ( gasaatr0(I).gt. gasco2) gasco2= gasco2tr0(I)
         if ( gasHNO3tr0(I).gt. gasHNO3) gasHNO3= gasHNO3tr0(I)
         if ( gasHCLtr0(I).gt. gasHcl) gasHcl= gasHcltr0(I)
         if ( gaslactr0(I).gt. gasLA) gasLA= gasLActr0(I)
         if ( gasoatr0(I).gt. gasOA) gasOA= gasOAtr0(I)
      enddo
c         read(1,*,end=222) timetr0(I), rhtr0(I),ttr0(I),
c     &  gasaatr0(I),gasamtr0(I),gasco2tr0(I),gashno3tr0(I),gashCLtr0(I)
c     & ,gasOAtr0(I) ,gaslactr0(I)
      if ( ML(3)+ gasaa.gt.0) then
         istrue(3)=1
         istrue(8)=1
         istrue(9)=1
      endif

      if ( ML(4)+ gasam.gt.0) then
         istrue(4)=1
         istrue(12)=1
         istrue(23)=1
      endif

      if ( ML(5)+ gasco2.gt.0) then
         istrue(5)=1
         istrue(13)=1
         istrue(14)=1
         istrue(15)=1
      endif
      xm24 =ml(24)+ml(25)+ml(26) ! phosphate
      if (xm24.gt.0d0) then
         istrue(24)=1
         istrue(25)=1
         istrue(26)=1
      endif
      xm27 =ml(27)+ml(28)
      if (xm27.gt.0d0) then
         istrue(27)=1
         istrue(28)=1
      endif
      
c     HNO3
      if ( ML(18)+ gasHNO3.gt.0d0) istrue(18)=1
c     hcl
      if ( ML(17)+ gasHcl.gt.0d0) istrue(17)=1
c     OA
      xm29=ml(29) +ml(30)+ml(31)
      
      if (xm29+gasoa.gt.0d0) then
         istrue(29)=1
         istrue(30)=1
         istrue(31)=1
      endif

      xm33=ml(33) +ml(34)
      if (xm33+gasLA.gt.0d0) then
         istrue(33)=1
         istrue(33)=1
      endif

      if (imode_ph.eq.2) then
         istrue(6)=0
         istrue(7)=0
      else
         istrue(6)=1
         istrue(7)=1
         
      endif
      
      DO j=1,np
         print*,j, ml(j), istrue(j)
      enddo


      

      iskinetic0=iskinetic
         iskinetic         =0
           IS=1
           NSS=1
           print*,'begin calH'
           
           call calHNew(T,mL)
c         DO I=1,NP
c         print*,i,ml(i)
c      enddo

           call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         
           print*, ' pH = ', dlog(ML(6)*gammah)/dlog(0.1d0),t
           print*,imode_NH4NO3,imode_MA
           print*, ml(6), gammah
           

           
           Do iI=1,np
           M(II)=ml(II)
           enddo
           DO kk=1 , 20
              ff=kk*.1
           Do iI=1,np
              ml(ii)=m(ii)*ff
              enddo
              ml(1)=m(1)
           call calHNew(T,mL)
           call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           print*, ml(29),aw, gamma2(1,29)
           
              enddo

           Do iI=1,np
              ML(II)=M(II)
           enddo
          
         iskinetic=iskinetic0

         call calHNew(T,mL)
         DO I=1,np
            print*,i,ml(I)
            
         enddo
         call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
       call vapnew(T,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
         
           print*, ' pH = ', dlog(ML(6)*gammah)/dlog(0.1d0)

           print*, 'aw', aw
           print*, 'CO2 ppm', PCO2 *1000
           

           
           print*,'g', gamma2(1,6), gamma2(1,7)
           print*,'g', gamma2(1,13), gamma2(1,14), gamma2(1,15)

c     calculate diffusion coefficient


              DO i=1, 200
                 ff=i/10d0
                              DO kk=1,np
                                 M(kk)=ff*ml(kk)
                                 enddo
              call aw_back
     & (t,M,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         call DH2O_NaCO3_cl(T,aw,m,dl,dlion)
c         write(52,*) aw, dl,dlion
c         write(54,*)  aw,dlion
c        if (aw.ge..5)  write(53,*) aw,dlion/dl

      enddo

      


      
         if (imode_output.eq.4) ML(21)=1D-6
         if (imode_output.eq.4) ML(22)=1D-6
         if (imode_output.eq.3) ML(21)=1D-6
         if (imode_output.eq.3) ML(22)=1D-6
         

         close(1)
         DO kk=1,np
            M(kk)=ML(kk)
            enddo

c           print*, Aw, poa
c         ml=0
 

c     adjust pH for mode 4  ( Adjust X+ and/ot Y-) to achieve the required pH at given gas phase, e.g. pH 6.6 of 2.5% CO2 and 135 ppb NH3)
          
         if (imode_output .eq. 4) then
c     calculate only one shell
            NSmax=1
            read(55,*) d19
            read(55,*) d20
            ML(19) = ml(19)*(1+d19)
            ML(20) = ml(20)*(1+d20)
            endif

            NSMAX0=NSMAX
            
 

c     initializing xn: th e moles of each species in each shell

c     set equal volume in each shell
c         vv=4*pi/3*x(2)**3
            time=0d0
c            rcore= core(time)
            rliq=r0
            
            call cal_rcore(time,rliq,rcore)
            
            call cal_MV(RH,t,mv)
            
            call set_shells_vol(T,RH,r0,ML,x,NS,nsmax)
c     recalculate xn21 and xn22
c     
            print*, rcore
            

            
         NSS=NS
         

               DO I=1,NS
                  ML(1)=1000d0/MM(1)
                  DO J=2,NP
                     ML(J)=ML(1)*xn(I,j)/xn(I,1)


                  enddo

                  IS=I
c                      print*,'is1 =' ,is
       call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        print*, 'aw,ph',I, aw, dlog(ml(6)*gammah)/dlog(.1d0)
        print*, 'm6', ml(6), gammah
        
               enddo
         
         

               
         nss=NS
         

c     set the boundary for mL(6) to a limited range of the previous value
c     in order to save CPU time

            is=1


              IS=NS
              
c                      print*,'is1 =' ,is
              DO J=2,NP
                 ML(J)=xn(NS,j)*ML(1)/xn(NS,1)
                 print*,J, ML(j)
              enddo

              call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        


       call vapnew(T,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
       print*, 'CO2 ppm ', pco2 *1000d0
                print*,'ph', NS, aw, dlog(ml(6)*gammah)/dlog(.1d0)
c                print*,'g', gammah, gamma2(1,6)
                

         gaa=0d0
         goA=0d0
         gLA=0d0
         ghcl=0d0
         gHNO3=0d0
         gco2=0d0
         gnh3=0d0
c     gas phase
         DO I=1,Ntr0
            if (gasaatr0(I).gt.gaa) gaa=gasaatr0(I)
            if (gasamtr0(I).gt.gNH3) gNH3=gasamtr0(I)
            if (gasco2tr0(I).gt.gco2) gco2=gasco2tr0(I)
            if (gashno3tr0(I).gt.ghno3) ghno3=gashno3tr0(I)
            if (gasoatr0(I).gt.goa) goa=gasOAtr0(I)
            if (gasLactr0(I).gt.gLA) gLA=gasLActr0(I)
            if (gasHCLtr0(I).gt.ghcl) gHCL=gasHCLtr0(I)
         enddo
         
c     condensed phase
            DO J=1,NP
               xn1(J)=0d0
            enddo

            DO J=1, NP
               DO I=1,ns
                  xn1(J)=xn1(J)+ xn(I,j)
               enddo
            enddo

c     add the gas phase

            
            xn1(3)=xn1(3)+ gaa
            xn1(4)=xn1(4)+ gNH3
            xn1(5)=xn1(5)+ gCO2
            xn1(17)=xn1(17)+ gHCL
            xn1(18)=xn1(18)+ gHNO3
            xn1(29)=xn1(29)+ gOA+xn1(30)+xn1(31)
            xn1(33)=xn1(33)+ glA+xn1(34)

            DO J=1, NP
               if (xn1(j).le.0d0) iszero(J)=1
               if (xn1(j).gt.0d0) iszero(J)=0
               if (dl_factor2(1,j).lt.0d0) iszero(J)=1
               
            enddo
c
           
            iszero(6)=1
            iszero(7)=1
            if (iskinetic0.eq.1) then
               iszero(13)=1
               iszero(5)=1
               iszero(15)=1
            endif

            
            if (iskinetic0.eq.0) then
               iszero(13)=1
               iszero(14)=1
               iszero(15)=1
            endif
            iszero(22)=1
            iszero(21)=1

            iszero(8)=1
            iszero(9)=1
            iszero(10)=1
            iszero(11)=1 ! protein big organics
            iszero(12)=1
            
            iszero(25)=1
            iszero(26)=1

            iszero(28)=1
            
            iszero(30)=1
            iszero(31)=1

            iszero(34)=1

            if (gasHNO3tr0(ntr0).gt.0d0) iszero(18)=0
            if (gasamtr0(ntr0).gt.0d0) iszero(4)=0
            if (gasaatr0(ntr0).gt.0d0) iszero(3)=0
            if (gasco2tr0(ntr0).gt.0d0) iszero(5)=0
            if (gasOAtr0(ntr0).gt.0d0) iszero(29)=0
            if (gaslactr0(ntr0).gt.0d0) iszero(33)=0
            if (gasHCLtr0(ntr0).gt.0d0) iszero(17)=0
            
            
            DO J=1, NP
               print*,'iszero',j, iszero(j),xn1(J)
            enddo

            
         DO I=1,NS
            xn(I,21)= 1d0/NS*6.023D-23 
            xn(I,22)= 1d0/NS*6.023D-23 
            enddo
c     total 1 virus 

         stiter0=0
         DO J=1,Ns
            stiter0=stiter0+xn(j,21)
          enddo
            
         vv=4*pi/3*x(2)**3
         ctiter0=xn(1,21)/vv ! initial virus concentration


      partvap3=gasaatr0(1)*press0*1D-9 ! hPa
      partvap4=gasamtr0(1)*press0*1D-9 ! hPa
      partvapco2=gasco2tr0(1)*press0*1D-6 ! hPa
      partvapHNO3=gasHNO3tr0(1)*press0*1D-9 ! hPa
      partvapHCL=gasHcltr0(1)*press0*1D-9 ! hPa
      partvapOA=gasOAtr0(1)*press0*1D-9 ! hPa
      partvaplac=gaslactr0(1)*press0*1D-9 ! hPa
      
c     
      partvapHNO3cond = 0
      partvap4cond = 0


      wtaero = 0
      
      
      DO I=1,ns
      wtaero = wtaero+      mv(2)* xn(I,2)
         DO J=7,np
            wtaero = wtaero+ mv(j)*xn(I,j)
         enddo
      enddo
c     calculate the H2O at ambient RH
      

      wtaero = wtaero * 1D6     ! mirogram / particle

      rh=.25                    ! RH of SMPS
      
      Do iI=1,np
      M(ii)=ML(ii)
      enddo
      DO jj=1,20000
         ff=jj* 0.02d0
         Do iI=1,np
         
         ml(II)=ff*M(II)
         enddo
         ml(1)=m(1)
c       call vapnew(T,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
c       if (a
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        if (aw.le.rh) goto 444
      enddo
      
 444  print*,'rh SMPS', rh,aw,wtaero

      DO I=1,NS
c     add volume of H2O
        xn11 = (xn(I,2)+xn(I,16)+xn(I,6))/(ml(2)+ml(16)+ml(6))*ml(1)
         wtaero = wtaero+      mv(1)* xn11*1E6
      enddo
c      print*,aw,wtaero

      
            xn_aero = pm/wtaero ! particles/m3

            xnHNO3c=0
            xnnh3c=0
            xnhclc=0
            xnp=0
            xnaac=0
            xnoac=0
            xnlacc=0
            print*, 'xn_aero ' , xn_aero
            DO I=1,NS
            xnNH3c =xnnh3c+
     & (xn(I,4)+xn(I,np-npsolid+1)+2*xn(I,np-npsolid+1))* xn_aero !moles /m3

            xnHNO3c = xnHNO3c+xn(I,18)* xn_aero !moles /m3
         xnHclc = xnHClc+(xn(I,17)+xn(I,npliq+9)+xn(I,npliq+7))* xn_aero !moles /m3 NaCl KCl
            xnp=xnp+(xn(I,24)+ xn(I,25)+xn(i,26)+xn(i,37) )*xn_aero
            xnaac=xnaac+xn(I,3)*xn_aero
            xnoac=xnoac+xn(I,29)+xn(I,30)+xn(I,31)
            xnlacc=xnlacc+xn(I,33)+xn(I,34)
            enddo

c     sum oxalic aicd in solids
            DO       I=1,NS
               DO Ip=1, 5
                        xnoac=xnoac + xn(I, ip+npliq)
               enddo
            enddo
            
            
            xnoac=xnoac*xn_aero

      partvapHClcond =  xnHCLc*8.314*T /100d0 ! hPa
      partvapHNO3cond =  xnHNO3c*8.314*T /100d0 ! hPa
      partvap4cond =  xnNH3c*8.314*T /100d0 ! hPa
      partvap3cond =  xnAAc*8.314*T /100d0 ! hPa
      partvapoacond =  xnoac*8.314*T /100d0 ! hPa
      partvaplaccond =  xnlacc*8.314*T /100d0 ! hPa

c      ppbP =  xnp*8.314*T /100d0/press*1D9 ! hPa
      ppbcl =  xnHClc*8.314*T /100d0/press*1D9 ! hPa
      print*, ' phosphor  ppb', ppbP
      print*, ' ppb cl ', ppbcl
      print*, ' condens, HNO3 ppb',partvapHNO3cond*1E6
      print*, ' condens, NH3 ppb',partvap4cond*1D6

      partvap4tot = partvap4 +partvap4cond
      partvap3tot = partvap3 +partvap3cond
      partvapHNO3tot = partvapHNO3 +partvapHNO3cond
      partvapHCltot = partvapHCl +partvapHClcond
      partvapoatot = partvapoA +partvapoacond
      partvaplactot = partvaplac +partvaplaccond

      
              IS=NS
              
c                      print*,'is1 =' ,is
              DO J=2,NP
                 ML(J)=xn(NS,j)*ML(1)/xn(NS,1)
                 print*,J, ML(j)
              enddo
        call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        

        print*, NS, aw, dlog(ml(6)*gammah)/dlog(.1d0)



       time=timetr0(1)
       print*, ' ns  =  ', ns,time



              NOUT=1
           print*, ' The first output may take several minutes to 1 h'

          INQUIRE (FILE='output_virus.dat', exist=ex)
          print*, 'ex', ex
              open(33,file='output_virus.dat')
              DO I=1,2
                 read(33,*,end=334,err=334) dd
              enddo
 334          N=I-1
              if (N.eq.0) ex=.FALSE.
              close(33)
              
              open(36,file='output_NH4NO3.dat')
              
c     Here read the data for restart  
c     open output file
              open(3,file='output_shells_M.dat')

c     Moles of all species in each shell
              open(23,file='output_shells_n.dat')
              open(21,file='output_shells_Dl.dat')
              open(34,file='output_n.dat')

c     active number of virus etc
              open(33,file='output_virus.dat')
c     partial pressures
              open(15,file='output_partial.dat')

              open(14,file='output_vapour.dat')
c     Molalities of all species in each shell

              open(18,file='output_area.dat')
           
              

c     continue the previous calculation if output_virus.dat exists
c     for a recalculation from the beginning, delete output_virus.dat
              IS=NS
              
c                      print*,'is1 =' ,is
              DO J=2,NP
                 ML(J)=xn(NS,j)*ML(1)/xn(NS,1)
                 print*,J, ML(j)
              enddo
        call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        print*,'befor ex'
        print*, NS, aw, dlog(ml(6)*gammah)/dlog(.1d0)


        
        if (ex ) then
      else
         if (isclose.ge.1) then
               open(35,file='output_total.dat')
      write(35,'(14E15.6)') 0d0, partvap4tot,partvapHNO3tot,partvapOAtot
     &        , partvapHCltot,partvap3tot, partvaplactot
      close(35)
      endif
      endif



                 if (ex ) then
c     set to the input value

c            DO I=1,NP-NPsolid
c                  print*,i,ml(I)
c               enddo
            
            
            NSMAX=NSMAX0
            print*,ns,nsmax
c
            
            
            if (isclose.ge.1) then
               open(35,file='output_total.dat')
               read(35,*)  xxx, partvap4tot,partvapHNO3tot,partvapOAtot
     &        , partvapHCltot,partvap3tot, partvaplactot
               close(35)
            endif

                 DO Ii=1,5000
                    goto 17
 18                 continue
                    print*,' conversion err'
 17                  read(33,*,end=19,err=19)   
     &     time,sflu,ssars,xnn,t
                    
                  read(14,*,end=19)timee
                  read(36,*,end=19)timee
                  read(15,*,end=19)timee,rh
                  read(34,*,end=19)timee

                  if (rh.le.0.27 ) ilate=1
                  
                  
                  enddo

 19                continue
                   print*,'inc0', inuc
                   DO I=1,Nnuc
                      
                      if (timee.ge.timeeff(I)) then
                         Inuc =inuc+1
                         print*,timee,timeeff(I)
                      endif
                      enddo
                      print*, 'inuc',timee, inuc,timeeff(inuc)

                   
                      if (ssars.le.1.1D-99 .and.
     &  imode_EDB+imode_NH3NO3+imode_MA.le.0  ) then
         print*,' finito, all the virus were inactivated stop'
         
         stop
      endif
      close(33)
      open(33,file='output_virus.dat', access='APPEND')
       if (time.ge.timetr0(Ntr0)) then
         print*,' finito, endtime reached stop'
         stop
      endif
      
c       if (imode_output.eq.3)    read(25,*) TIMEe, partvap4tot
c     & ,  partvaphno3tot

       if (imode_output.eq.3) print*, 'total NH3', partvap4tot
       if (imode_output.eq.3) print*, 'total HNO3', partvapHNO3tot

                   Nout=Ii-1
                   if (nout.lt.1) Nout=1
                   
                   print*,time, outputtime(nout), outputtime(nout+1)
                  
                   print*, 'NOUT ', nOut
               DO Ii=1,NOUT
              read(3,*,end=155) time,x(NS+1),ns,x(1)

              DO I = 1,ns
              read(3,*,end=155) jj,
     &        x(I+1), awshell(i)

              enddo

              
           if (time.ge. .1*(r0*1D4)**2  .and.NS.ge.2.and.istimeeq.eq.0
     & .and. imode_output.ne.2                      ) then
         dawmax =0
        xkelvin = dexp( 2* sigma * MV(1) /(8.314E7*T*x(NS+1)) )

         DO Iii=1,ns
            dd=dabs(awshell(iii)- rhtr0(Ntr0)/xkelvin)
            if (dd.ge. dawmax) dawmax=dd
         enddo
         if (dawmax.le.0.01d0) then
            timeeq= time
            istimeeq=1
         endif
         endif
    
c         print*, timeeq, 'timeeq',ns
         
      


              DO I = 1,ns
              read(23,*,end=155) Iii,time11, (xn(i,j),j=1,NP)
              read(21,*,end=155) Iii
              enddo

              
            enddo
              goto 156
 155          Nout=II-1
              print*,'NOUT =', Nout

 156          continue

              DO I = 1,Noutput
                  if (outputtime(I).gt. time11) goto 157
                  enddo
 157              Nout=I
                  print*,NS, timeeq

                  
c     if iscenter =1 set all solid into shell 1
c     it useful when calculate for very long time, set imode_output > 10 and reduce the number of shells to make rapid runs.
              DO I=1,NS
              IS=I
              
              DO J=2,NP
                 mL(J)= ml(1)*xn(i,J)/xn(i,1)
c                 print*, j,ml(J)
              enddo
              
c           call calHNew(T,mL)
           call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
c
c         cccccccccccccccccccccccc
c     calculate the activity coefficients for the first time
           
       gamma2(I,13) =get_gammaco2(t,Ml)       
       gamma2(I,6)=gammah
       gamma2(I,12)=gammaNH4
       gamma2(I,17)=gammaCl
       gamma2(I,16)=gammaNa
       gamma2(I,18)=gammaNO3
       gamma2(I,27)=gammas1
       gamma2(I,28)=gammaS2

       gamma2(I,1)=aw
       xmplus=ml(6)+ml(19)+ml(12)+ml(16)
       gammaOH =get_gammaOH(xMplus)
       gamma2(I,7)=gammaoH

       gamma2(I,15)=gammaHCO3
       gamma2(I,14)=gammaCO3

       gamma2(I,24)=gammaH3PO4
       gamma2(I,25)=gammaH2PO4
       gamma2(I,26)=gammaHPO4

       gamma2(I,29)=gammaH2OA
       gamma2(I,30)=gammaHOA
       gamma2(I,31)=gammaOA

       gamma2(I,19)=gammak
       gamma2(I,36)=gammaca
       gamma2(I,32)=gammamg
cccccccccccccccccccccccc
       print*,'aw,pH', i,aw, dlog(gammah*ml(6))/dlog(.1d0)

      enddo

              

              if (iscenter.eq.1 ) then
              DO I=2,NS
                 DO J=npl+1, NP
                 xn(1,j)= xn(1,j)+ xn(I,j)
                  xn(I,j)=0d0
                    enddo
                    enddo
                    endif


c     For restart calculation, the number of shells can be less or equal than previus runs, but not more.
c     
              
c     for restart, one can reduce the number of shells, but not increase!
c     It is useful to reduce the number of shells when equilibrium is reached and the calculation is slow due to high number of shells and one want to reach longer times (months - years)
              if (imode_shell.eq.1)nsmax=ns
                print*,ns,nsmax
                xn22=0
                DO I=1,NS
                   xn22=xn22+ xn(I,2)
                enddo
                print*,'xn2', xn22

                
                
c     group to 1 bin
              if(NSmax.eq.1 .and.NS.gt.1) then

                 DO I=2, NS
                    DO J=1,NP
                       xn(1,j) =xn(1,j) + xn(I,j)
                    enddo
                    enddo

                 NS=1
                 endif
     
c     reduces shells  to NSmax shells
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &        +xn(I,npl+11)      +xn(I,37)
         

         enddo

         print*, 'xnp, befor regroup', xnp,nsmax,ns
                 
         
                 if (NSmax.le.NS/2) then
c     

                    DO I=1,NSmax
                          DO J=1,NP
                             xnnew(i,j)=0d0
                          enddo
                          enddo

                    N = NS/NSMAX
                    Nrest= NS- N*NSMAx

                    N2=0
                    
                    DO I= 1, NSMAx
                       NN=N
                       if (I.le.Nrest)  NN=N+1
                       N1= N2+1
                       N2= N2+NN
                       print*,'n1,n2', N1,N2,ns
                       DO II=N1,N2
                          DO J=1,NP
                             xnnew(I,j) =xnnew(I,j)+ xn(Ii,j) 
                          enddo
                       enddo
                    enddo


                    
                    NS=NSMAx

                    DO I=1,NSmax
                          DO J=1,NP
                             xn(I,J)=xnnew(i,j)
                          enddo
                          print*, I, xn(I,1),xn(I,16)
                          enddo

                       endif
                                             

c     group outer  bins into one bin, when  2/NS < Nsmax < NS
                    print*,'nn',ns/nsmax
                    
                    if(ns/nsmax.eq.1 .and. NSmax.lt.ns) then
                       print*,'nsmax > ns ', nsmax,ns
                       

                       DO I=NSmax+1, Ns
                         
                             DO J=1,NP
                                xn(nsmax,J) =xn(nsmax,J) +xn(I,J) 
                                enddo
                             enddo

                          NS=NSMAx
                       endif
                       

c       if NSmax is a multiple of NS
                       NN=nsmax/NS
                       if (nn* NS.eq. nsmax .and. nn.gt.1)then
                          
                          DO I=1,ns
                             DO II=1,NN
                                DO J=1,Np-npsolid
                                   xnnew(nn*(I-1)+ii,j) = xn(I,j)/nn
                                enddo
                                DO J=Np-npsolid+1,np

                                   xnnew(nn*(I-1) +ii,j) = 0d0
                           if (II.eq.1) xnnew(nn*(I-1) +ii,j) = xn(I,j)

                         enddo
                             enddo
                       enddo
                      NS=NSmax
                      DO I=1,NS
                         DO J=1,np
                            xn(i,j)= xnnew(i,j)
                         enddo
                         enddo
                         
                       endif
                       

                          print*,' ns= ', ns
               NOUT=NOUT+1
               
               DO I = 1,Noutput
                  if (outputtime(I).gt. time+1D-6) goto 119
               enddo
 119           print*, 'time', time, outputtime(I)
               nout=I
               print*,'ns',ns
               
               DO I=1,NS
                  ML(1)=1000d0/MM(1)
                  DO J=2,NP
                     ML(J)=ML(1)*xn(I,j)/xn(I,1)

c                  print*,j,ml(J)

                  enddo
                      IS=I
c                      print*,'is1 =' ,is
       call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        print*, I, aw
               enddo


               I=Ns
               print*,(x(j),j=1,ns+1)
               print*,t,time
               print*, ml(4),ml(29)
               

               

               



       call vapnew(T,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
       print*,'aw',aw,pnh3
       
       
c       print*, 's_pp ', pnh3*poa*dsqrt(aw)/apNH4HC2O4_05h2o(T)


       
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo

         
         
             endif
             
cccccccccccccccccccccccccccccccccccccccccccccccc
c     start of the time integral
cccccccccccccccccccccccccccccccccccccccccccccccc
c     test

               DO I=1,NS
                  ML(1)=1000d0/MM(1)
                  DO J=2,NP
                     ML(J)=ML(1)*xn(I,j)/xn(I,1)

c                  print*,j,ml(J)

                  enddo
                      IS=I
c                      print*,'is1 =' ,is
                   iskinetic=0
                      
       call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        print*, I, aw, dlog(gammah*ml(6))/dlog(.1d0)

        iskinetic=iskinetic0

        call calHNew(T,mL)
        call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        print*, I, aw, dlog(gammah*ml(6))/dlog(.1d0)

        print*, ml(5),ml(13),t

c        xn(I,13) = mL(1)/xn(I,1)* ml(13)
c        xn(I,14) = mL(1)/xn(I,1)* ml(14)
c        xn(I,15) = mL(1)/xn(I,1)* ml(15)
        

      enddo


      
      IS=NS
      DO jj=1,np
      M(jj)=ml(jj)
      enddo
        call cal_ml(T,rh,ML)
      print*, 'm18',m(18),ml(18),rh
      
      call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        print*, 'after cal_ML', aw, dlog(gammah*ml(6))/dlog(.1d0)

      
      dtimeb=1D-10
      dtime=dtimeb
      
      print*, 'begin time loop'
      nmin=ns
      jmin=1001
            call flush(6)
c            istimeeq=0

             xn22=0d0
             DO I=1,ns
                xn22=xn22+xn(I,2)
             enddo
             print*, 'xn22 befor time loop', ns,xn22,imode_MA
             
             if (imode_MA.eq.1) then
                close(22)
c                DO I=1,ntr0
                   
                partvapoa=0d0
c                time=timetr0(I)
                pp= ppvap0 + (ppvap1-ppvap0)*(1-
     &            DEXP(-TAUP *
     &             (time-timetr0(1))/(timetr0(ntr0)-timetr0(1))
     &             ))
c                RHH=RHTR0(i)
                RHH=RH
                if (rhh.le.0.15) rhh=0.15

             pp=pp* (1+ dppvap/100d0*dlog(rhh/.4d0))
             partvapoa=pp*press0*1D-9 ! hPa
c             write(22,*) time ,pp
             
c          enddo
             
          endif
c             stop
             
         if (imode_NH4NO3.eq.1) then
             ppvap= ppvap0 + (ppvap1-ppvap0)*(1-
     &            DEXP(-TAUP *
     &             (time-timetr0(1))/(timetr0(ntr0)-timetr0(1))
     &             ))
             rhh=rh
             
             if (rhh.le.0.15) rhh=0.15
             ppvap=ppvap* (1+ dppvap/100d0*dlog(rhh/.4d0))
             
        call cal_flux(time,x,flux,dtime,NS)
        
      endif

c     get the final size
         M=ML
         M=0
         m(1)= ml(1)
         m(2)= ml(2)*2
         print*,'imode_NH4NO3', imode_NH4NO3
c         stop
         if (imode_NH4NO3.eq.1 .or. imode_MA) then
         open(19,file='output_drysize.dat')
         DO I=1,ntr0
         rh1= rhtr0(I)
           call cal_ml(T,rh1,M)

           call cal_MV(rh1,t,mv)

           xn11= xn22/m(2)*m(1)
         vol= xn22*mv(2)+ xn11*mv(1)
         rr= (vol/4/pi*3d0) **(1/3d0)

         write(19,*) timetr0(I), rr
      enddo
      print*,'finished drysize'
      endif
         close(19)
        
         dtime00=dtime
         
 11   continue

      dtime=dtime00
       if (time+dtime.gt. outputtime(Nout)) then
      dtime= outputtime(NoUT)-time
      if (dtime.le.1D-12)dtime=1D-12
      endif

      
c     rcore=core(time)
c           print*, time, rcore,x(1)
      rliq=(x(NS+1)**3-x(1)**3)**(1/3d0)

      call cal_rcore(time,rliq,rcore)

      x(1)=rcore
           
           if (isNAN(awshell(1))) then
           print*, 'NAN,dtime', dtime
           stop
        endif
c     print*, timeeq, 'timeeq'
c           print*,'dtimeee', dtime
           if (dtime.le.1D-19)then
              print*, 'dtime < 1E-19s stop',dtime
              stop
              
           endif
           if (idiagnose.ge.4) print*,'aaa9'

      if (time.ge. .01*(r0*1D4)**2 .and.istimeeq.eq.0 .and.NS.ge.2) then
         dawmax =0
        xkelvin = dexp( 2* sigma * MV(1) /(8.314E7*T*x(NS+1)) )

         DO I=1,ns
            dd=dabs(awshell(i)- rhtr0(Ntr0)/xkelvin)
            if (dd.ge. dawmax) dawmax=dd
         enddo
         if (dawmax.le.0.01d0) then
            timeeq= time
            istimeeq=1
         endif
         endif

c         print*,' timeeq',time, timeeq
         
c     no H2O equilibrium for output mode 2 (EDB)
         if (imode_output.eq.2)  then
            timeeq= 1D20
            istimeeq=0
         endif
         
              
         rhh1=rh-.3
         if (rhh1.le.0d0) rhh1=0
         dxmin=dxmin0 * (1+ RHh1/.1)
         
              n1=1
          time1(1)=time
          timeA=time

c     define the partial pressure
          
c     interpolates the T

c     if Eddy Diffusion coefficient is zero, then make linear interpolation of input data
c
c          make precise pH calculation at the end
           
          if (imode_output.eq.3 .and. timetr0(NTR0)-time.le.100) 
     & imode_PH=0

          

           

          if (Deddy.le.1D-10 .or. time.lt.0d0) then
c     T
                  time1(1)=time
             call intpl(timetr0, ttr0, ntr0, time1, y1, n1)          
          
             
             t=y1(1)
c     RH
         call intpl(timetr0, gash2otr0, ntr0, time1, y1, n1)          
         rh=y1(1)/vwater(T)
c
         if(time.ge.0) then
         call intpl(timetr0, rhtr0, ntr0, time1, y1, n1)          
         rh=y1(1)

         endif

         if (imode_NH4NO3.eq.1) then
             ppvap= ppvap0 + (ppvap1-ppvap0)*(1-
     &            DEXP(-TAUP *
     &             (time-timetr0(1))/(timetr0(ntr0)-timetr0(1))
     &             ))
             rhh=rh
             if (rhh.le.0.15) rhh=0.15
             ppvap=ppvap* (1+ dppvap/100d0*dlog(rhh/.4d0))
          endif

          

        call intpl(timetr0, gasaatr0, ntr0, time1, y1, n1)          
         ppbace=y1(1)
         call intpl(timetr0, gasamtr0, ntr0, time1, y1, n1)          
                  ppbnh3=y1(1)
               call intpl(timetr0, gasco2tr0, ntr0, time1, y1, n1)          
                  partvapco2=y1(1)*1D-6*press0
                  call intpl(timetr0, gashno3tr0, ntr0, time1, y1, n1)          
                  partvaphno3=y1(1)*1D-9*press0  
                  call intpl(timetr0, gashcltr0, ntr0, time1, y1, n1)          
                  partvaphcl=y1(1)*1D-9*press0   
                  call intpl(timetr0, gasOAtr0, ntr0, time1, y1, n1)          
                  partvapOA=y1(1)*1D-9*press0   

                  if (imode_MA.eq.1) then
                pp= ppvap0 + (ppvap1-ppvap0)*(1-
     &            DEXP(-TAUP *
     &             (time-timetr0(1))/(timetr0(ntr0)-timetr0(1))
     &             ))
             rhh=rh

             if (rhh.le.0.15) rhh=0.15
             pp=pp* (1+ dppvap/100d0*dlog(rhh/.4d0))
             partvapOA=pp*press0*1D-9 ! hPa

          endif
         

                  call intpl(timetr0, gaslactr0, ntr0, time1, y1, n1)          
                  partvaplac=y1(1)*1D-9*press0   

                  
               else
c     take the eddy diffustion coefficient of time at 0 and the 
                  
                  A0= 2*pi*0.75**2
                  At = a0+ 4 * Deddy *time
                  if (time.le.0) At= A0
                  
            t = a0/At * ttr0(iexit) + (At-a0)/At*ttr0(Ntr0)
              pwater = a0/At * gash2otr0(iexit)
     & + (At-a0)/At*gash2otr0(ntr0)
              RH=pwater/vwater(T)

         ppbace = a0/At * gasaatr0(iexit) + (At-a0)/At*gasaatr0(Ntr0)
         ppbNH3 = a0/At * gasamtr0(iexit) + (At-a0)/At*gasamtr0(Ntr0)
         xxx = a0/At * gasco2tr0(iexit) + (At-a0)/At*gasco2tr0(Ntr0)
         partvapco2= xxx*1D-6*press0


         xxx = a0/At * gasHNO3tr0(iexit) + (At-a0)/At*gasHNO3tr0(Ntr0)
         partvaphno3= xxx*1D-9*press0         

         xxx = a0/At * gasHcltr0(iexit) + (At-a0)/At*gasHcltr0(Ntr0)
         partvaphcl= xxx*1D-9*press0         
         xxx = a0/At * gasOAtr0(iexit) + (At-a0)/At*gasOAtr0(Ntr0)
         partvapOA= xxx*1D-9*press0         

         xxx = a0/At * gaslactr0(iexit) + (At-a0)/At*gaslactr0(Ntr0)
         partvaplac= xxx*1D-9*press0         

         endif


       if (rh.le. .01d0) rh=.01d0 !not dryer than 1 %
         partvap3=ppbace*1D-9*press0
         partvap4=ppbnh3*1D-9*press0
         
          call intpl(timetr0, dgasco2tr0, ntr0, time1, y1, n1)          
          isdegas=1
          if (y1(1).gt. 1D-20) isdegas=0d0
          

c     apply only to imode_output 3
c     calculate particle pressure for NH3  or HNO3 when it is closed system

c     by apply HNO3 od oxaliy acid into gas phase, NH3 goes to the background aerosol. NH3 is treated as a closed system
c     
c     by removal NH3, HNO3 may release from aerosol


          
          if (isNH3close+  isHNO3close + isOAclose
     & +isAAclose+ ishclclose+islaclose.ge.1 ) then

                xnNH3c=0
                xnHNO3c=0
                xnOAc=0
                xnAAc=0
                xnHClc=0
                xnlacc=0
                
                do i=1,ns
                   
               xnNH3c =XnNH3C+
     & ( xn(i,4)+2*xn(i,np-npsolid+2)+ xn(i,np-npsolid+1))* xn_aero !moles /m3
                   xnOAc =XnOAC+
     &(xn(i,29)+xn(i,30)+xn(i,31)+xn(i,np-npsolid+1)+xn(i,np-npsolid+2)+
     & xn(i,np-npsolid+3)+xn(i,np-npsolid+4)+xn(i,np-npsolid+5))*xn_aero !moles /m3
         xnHclc=xnHclc+(xn(I,17)+xn(i,npl+9)+xn(i,np-npsolid+7))*xn_aero  !moles /m3
            xnAAc = xnAAc+(xn(I,3))* xn_aero  !moles /m3
            xnlacc = xnLAcc+(xn(I,33)+xn(I,34))* xn_aero  !moles /m3
            xnHNO3c = xnHNO3c+xn(I,18)* xn_aero  !moles /m3
         ENDDO
         
            partvapHNO3cond = xnHNO3c*8.314*T /100d0 ! hPa
            partvapHclcond =  xnHclc*8.314*T /100d0 ! hPa
            partvap4cond =  xnNH3c*8.314*T /100d0 ! hPa
            partvapOAcond =  xnOAc*8.314*T /100d0 ! hPa
            partvap3cond =  xnAAc*8.314*T /100d0 ! hPa
            partvaplaccond =  xnlacc*8.314*T /100d0 ! hPa

            
            if (isNH3close.eq.1) then
               if (partvap4cond.gt.partvap4tot) then
                  ff= partvap4tot/partvap4cond
                  partvap4cond=partvap4tot
                  DO I=1,NS
                     xn(I,4) =ff*xn(I,4)
                  enddo
               endif
               
               partvap4 =partvap4tot -partvap4cond
            endif

            if (isAAclose.eq.1) then
               if (partvap3cond.gt.partvap3tot) then
                  ff= partvap3tot/partvap3cond
                  partvap3cond=partvap3tot
                  DO I=1,NS
                     xn(I,3) =ff*xn(I,3)
                  enddo
               endif
               
               partvap3 =partvap3tot -partvap3cond
            endif

            if (islacclose.eq.1) then
               if (partvaplaccond.gt.partvaplactot) then
                  ff= partvaplactot/partvaplaccond
                  partvaplaccond=partvaplactot
                  DO I=1,NS
                     xn(I,33) =ff*xn(I,33)
                     xn(I,34) =ff*xn(I,34)
                  enddo
               endif
               
               partvap3 =partvap3tot -partvap3cond
            endif
            
            
            if(isHNO3close.eq.1) then
               if (partvapHNO3cond.gt.partvapHNO3tot) then
                  ff= partvapHNO3tot/partvapHNO3cond
                  partvapHNO3cond=partvapHNO3tot
                  DO I=1,NS
                     xn(I,18) =ff*xn(I,18)
                  enddo
               endif
               

               partvapHNO3=partvapHNO3tot-partvapHNO3cond
               endif

            if(isHClclose.eq.1) then
               if (partvapHClcond.gt.partvapHCltot) then
                  ff= partvapHCltot/partvapHClcond
                  partvapHClcond=partvapHCltot
                  DO I=1,NS
                     xn(I,17) =ff*xn(I,17)
                     xn(I,npl+9) =ff*xn(I,npl+9)
                     xn(I,npl+7) =ff*xn(I,npl+7)
                  enddo
               endif
               

               partvapHcl=partvapHcltot-partvapHclcond
               endif

               if(isOAclose.eq.1)then
               if (partvapOAcond.gt.partvapOAtot) then
                  ff= partvapOAtot/partvapOAcond
                  partvapOAcond=partvapOAtot
                  DO I=1,NS
                     xn(I,29) =ff*xn(I,29)
                     xn(I,30) =ff*xn(I,30)
                     xn(I,31) =ff*xn(I,31)
                  enddo
               endif

                  partvapoa=partvapoatot-partvapoacond
                  endif
      endif


              press=press0

              xnsolidt=0
              DO kk=1,ns
              xnsolidshell(kk)=0
              enddo
              DO I=1,NS
                 xnsolidt=xn(I,np)+xnsolidt
                 DO kk=1,npsolid
                 xnsolidshell(i)=xn(I,kk+np-npsolid)+xnsolidshell(i)
                 xnsolidt=xn(I,kk+np-npsolid)+xnsolidt
                 enddo
              enddo
c     end define partial pressures
              xnsolid=xnsolidt

              
c     merge shell with less molecules to the neighboring shells
c         if ( xnsolid.le.1D-30) then
              
c           if (isph.eq.0) then  
c              call reset_shells(NS,x)
c           else
              
c           endif
           
c     else
c          call reset_shells_solid(NS,x)
c        endif
c        print*, '  pre cal'
            sigma= cal_sigma(Tdrop,RH)
c     for low pH < 3 now kinetic anymore
            sphh=0
            s5=0
            s13=0
            
            DO I=1,NS
               sphh=sphh+phshell(I)/NS
               s13=s13+ xn(I,13)
               s5=s5+ xn(I,5)
            enddo
c     now kinetic calculation when pH < 4 and few carbonate
            dliq=x(NS+1)-x(1)
            
            
            iskinetic=iskinetic0
c            if ( s5.le.1E-30) iskinetic=0
c            if ( s5.le.1E-30) timeeq=timeeq0

            if ( s5.ge.1E-30) then
c               if (sphh.le.4 .and. s13/s5.ge.0.5d0)  iskinetic=0
c               if (sphh.le.4 .and. s13/s5.ge.0.5d0)   timeeq=timeeq0
            endif
            


            if(Ns.eq.1 .and. time.ge.1) then
               xis=x(NS)
               vv= 4*pi/3*x(NS)**3
        DO kk=np-npsolid+1,np
           vv=vv+ mv(kk)*xn(1,kk)
           enddo
           xis=(vv/4/pi*3)**(1/3d0)

           dx = x(NS+1)-xis
c
c     thickness < 1 nm H2O equilibrium
 
           if (dx.le.1D-7) timeeq=time-1

           dliq=dx
           
        endif

c        print*,'dliq', dx,dliq
        
       tdrop1=ta
c
c        if (idiagnose.ge.4) print*, 'aa1'
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo

      xnca=0
      DO I=1,ns

         xnca=xnca+ xn(I,36)+xn(I,npl+10)+xn(I,npl+11)         

         enddo
         
c     end solid
                                       
         if (idiagnose.ge.4) print*, 'aa1,pre cal_flux',xnca

c     jmin=-1
c        nmin=-1
         NN=NS-3
         if (NN.le.0) NN=1
        call cal_flux(time,x,flux,dtime,NS)
c           print*,'dtime a', dtime
c        print*, 'ishelleq af cal_flux', ishelleq
        
        if (idiagnose.ge.1) then
           
           write(91,'(16E15.6)') time, 
     &  (awshell(I),i=NS,1,-1)
       write(92,'(16E15.6)') time, (flux(I+1,1),i=NS,1,-1)

       write(93,'(16E15.6)') time, (xn(I+1,1),i=NS,1,-1)

      endif
      
c     print*, 'after cal_f', flux(NS+1,1)
c        if (idiagnose.ge.4) print*, 'aa2'

c     write(6,'(10E15.6)') time, (awshell(I),i=1,ns)
        
c        write(6,'(10E15.6)') time, flux(4,ns+1), flux(ns+1,18)

c     for NH4No3
c     take the same flux for both ions NH4+, NO3-1
c     it means, NH4 NO3 diffuse always in pairs
        
        if (imode_NH4NO3.eq.1) then
c     avoid small deviation fo flux
           DO I=2,NS+1
              dd=flux(I,4)+flux(I,18)
              flux(I,4)= dd/2
              flux(I,18)= dd/2
              
           enddo
           endif
c        write(6,'(A,10E15.6)') 'aa',time, flux(ns+1,4),flux(ns+1,18)
             
         tdrop=ta


        loop=loop+1
        xloop=xloop+1d0

c     calculate the time step and fluxes

        kka=0
c        loopkka=        loopkka+1
        isdiv=0
        
        DO I=1,NS
           is=I
           DO j=1,Np
              xn0(i,j)=xn(i,j)
              gamma20(i,j)=gamma2(i,j)
           enddo

           
              if (dabs(Tdrop-tdrop1).ge.1D13) then
         DO J=2,npliq
            ml(J)= ml(1)*xn0(I,j)/xn0(I,1)
         enddo
         IS=I
         tdrop1=tdrop
         call calHNew(Tdrop,mL)
          call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
       phshell(i) = dlog(ml(6)*gamma2(I,6))/dlog(.1d0)
       phshell0(i) = dlog(ml(6)*gamma2(I,6))/dlog(.1d0)
        awshell0(i) =aw
        awshell(i) =aw
c        print*,'ttt', i, dlog(ml(6)*gamma2(I,6))/dlog(.1d0),aw
      endif

        phshell0(i) =phshell(I)
        awshell0(i) =awshell(I)
      
      enddo

      
           
 499    continue

        kka=kka+1
        loopkka=        loopkka+1

c     roll back of the compostion
           DO I=1,NS
           DO j=1,Np
              xn(i,j)=xn0(i,j)
           enddo
           enddo
           
c        endif
        

           
        if (kka.ge.100)then
          write(6,'(A,2I5,14E15.6)') 'dpH',kka,Imax,
     & dphmax,phhmax,phshell0(Imax)
     *     ,dtime,time,tdrop,tdrop1
          
      print*,'error: dpH problem',time
           stop
        endif

c        print*,'dtime b',dtime
            


c     calculates the titer of viruses
c     grout the dissociate spieces
        if (kka.le.1) then
        DO I=1,NS-1
           flux(I+1,3)=flux(I+1,8)+flux(I+1,9)+flux(I+1,10)
           
           FLUX(i+1,8)=0D0
           FLUX(i+1,9)=0D0
           FLUX(i+1,10)=0D0

           FLUX(i+1,6)=0D0
           FLUX(i+1,7)=0D0

           flux(I+1,4)=flux(I+1,12)+flux(I+1,23)+flux(I+1,10)

           FLUX(i+1,12)=0D0
           FLUX(i+1,23)=0D0
c           print*,'flux3',flux(I+1,3), xn(I,3tr
           
           
        flux(I+1,24)=flux(I+1,24)+flux(I+1,25)+flux(I+1,26)+
     &         flux(I+1,37)
       flux(I+1,26)=0d0
       flux(I+1,25)=0d0
       flux(I+1,37)=0d0


           flux(I+1,27)=flux(I+1,27)+flux(I+1,28)
          flux(I+1,28)=0d0
       

       flux(I+1,29)=flux(I+1,29)+flux(I+1,30)+flux(I+1,31)+flux(I+1,38)
       flux(I+1,30)=0d0
       flux(I+1,31)=0d0
       flux(I+1,38)=0d0

       flux(I+1,33)=flux(I+1,33)+flux(I+1,34)
       flux(I+1,34)=0d0

       if (iskinetic.eq.0) then

           flux(I+1,5)=flux(I+1,13)+flux(I+1,14)+flux(I+1,15)
           flux(I+1,13)=0d0
           flux(I+1,14)=0d0
           flux(I+1,15)=0d0
        endif

        if (iskinetic.eq.1) then

           flux(I+1,5)=0d0
           flux(I+1,13)=0d0
           flux(I+1,15)=0d0
           flux(1,5)=0d0
           
        endif
       
      enddo
      if (iseqAA.eq.1) then
        flux(NS+1,3) =flux(NS,3)
      endif
      
      if (iseqNH3.eq.1) then
        flux(NS+1,4) =flux(NS,4)
      endif

      endif

      
      
        DO I=1,NS

           xn(I,3) =xn(I,8)+xn(I,9)+xn(I,10)
           xn(I,4) =xn(I,12)+xn(I,23)+xn(I,10)
           
           xn(I,24) =xn(I,24)+xn(I,25)+xn(I,26)+xn(I,37)
           xn(I,25)=0d0
           xn(I,26)=0d0
           xn(I,37)=0d0


           xn(I,27) =xn(I,27)+xn(I,28)
           xn(I,28)=0d0
       
           xn(I,29) =xn(I,29)+xn(I,30)+xn(I,31)+xn(I,38)
           xn(I,30)=0d0
           xn(I,31)=0d0
           xn(I,38)=0d0

           xn(I,33) =xn(I,33)+xn(I,34)
           xn(I,34)=0d0

       if (iskinetic.eq.0) then
           xn(I,5) =xn(I,13)+xn(I,14)+xn(I,15)
           xn(I,13)=0d0
           xn(I,14)=0d0
           xn(I,15)=0d0

        endif
       if (iskinetic.eq.1) then
           xn(I,14)= xn(I,14) +xn(I,15) ! take the sum of CO3 adn HCO3 ions
           xn(I,15)=0d0
        endif
       
      enddo
          
          


c     ccheckf flux
          II1=1
          if (imode_shell.eq.1) II1=2
          DO II=1,II1
          DO J=1, Np-npsolid
             

             
             
              if (dl_factor2(1,j).gt.0 .and. iszero(j).eq.0) then
                DO I=1,ns
                 if (imode_NH4no3.eq.1.and.I.eq.ns) then
                    if (J.eq.18 .or.j.eq.4) then
                    if ( ml(1)*xn(ns,j)/xn(ns,1).le.1D-12) then
                       flux(I+1,j)=flux(I,j)
                    endif
                 endif
                 
                 endif

                   df=(flux(I+1,j)-flux(I,j))
                  
             if (iskinetic.eq.0) then
                if (J.eq.13) df=0
                if (J.eq.14) df=0
                if (J.eq.15) df=0
                
                if (J.eq.5 .and. I.eq.NS) df=0d0 !equilibrium of CO2 last shell
             endif
             

                if (iskinetic.eq.1) then
                if (J.eq.5  ) df=0d0
             endif


                if (df .lt.0d0) then
                   if (xn(I,j) .gt.1D-51) then
                      dtime1=dtime
                      
              if (xn(I,j)+ df*dtime.lt.0.and.isdiv.le.1 ) then

                 dtime=dabs(xn(I,j)/df)/2

                 if (Idiagnose.ge.1) then
                 write(6,'(A,3I5,5E15.6)')'bb liq',
     &                 NS,i,j,xn(i,j)/xn(I,1)*ml(1),dtime,dtime1,
     &                flux(I+1,j), flux(I,j)
                 phh=dlog( xn(I,6)*ml(1)/xn(I,1))/dlog(.1d0)
c                 write(6,'(A,16E15.6)')'f gma',
c     &     gamma2(i,6),  gamma2(i,7),  gamma2(i,13),  gamma2(i,15)
c     & ,gamma2(I,14),gamma2(3,6),gamma2(3,14),phh
                 
              endif
                 nmin=I
          jmin= j


          
       endif

              else

                 print*,'warning: xn(i,j) <0) liq', i,j, xn(i,j)
                 print*, 'aa',xn(i,j)

              endif
              

              endif
                    
                  

              enddo
              
           endif

        enddo
        enddo
c        print*,'dtime c',dtime

c     endcheck
           if (isdiv.eq.2) goto 599


c     calculates the new composition after the timestep dtime
c     for dissociation species, the liquid phase diffusion is set to the sum of the individual species.

          

       isolidmin=0
         
ccccccccccccc DO solid part ccccccc

         
         DO I=1,ns
          DO KS = 1, npsolid
             ISs= NP-NPsolid+KS

          if (xn(I,NP-npsolid+KS) .gt.0d0) then



             ff=(flux(I,np+KS))+(flux(I,IsS))
             xnss=xn(I,iss)+ ff*dtime
             
             if (xnss.lt.0d0) then

                dtime1=dtime
                dtime =dabs(xn(I,iss) /ff)
                Imin=I
                isolidmin=KS
                j=0
                
                if (Idiagnose.ge.1) then
                   write(6,'(A,2I5,5E15.6)')'bb solid I,ks', I,KS,
     &xn(I,iss)*ml(1)/xn(I,1),
     &               xnss     *ml(1)/xn(I,1),dtime,dtime1
                   
                   
                endif

             endif

                

c flux to solids
        

c     check maximal flux
                    DO K=1,Nsp(KS),-1
c     

                       j1=nsp_index(KS,K) !species index
                       j=j1
                       if (j1.eq.25 .or. j1.eq.26 .or.j1.eq.24)j=24
                       if (j1.eq.29 .or. j1.eq.30 .or.j1.eq.31)j=29
                       if (iskinetic.eq.0) then
                          if (j1.eq.14 .or.j1.eq.15) J=5
                       endif

                       xnu=xnsp_nv(Ks,K) !how many

c     flux from solid  to shell i
                         xnss=xn(I, j) - xnu*flux(I,np+ks)*dtime

                         if (xnss.lt.0d0) then
                            dtime1=dtime
c     3.3: It will be 10% left even 3 solids reduce the same species
c                            dtime=xn(I, j)/flux(I,np+ks)/xnu/3.3 !too big change of liquid
                 if (Idiagnose.ge.1) then
                      write(6,'(A,4I5,5E15.6)')'bb solid I',
     &      i,ks,j,j1,xn(i,j)/xn(I,1)*ml(1),dtime,dtime1

                            endif

          nmin=I
          jmin=  J
              isolidmin=ks


                       endif
                         
                           if (I.ge.2) then
c     flux from solid  to shell i-1
                         xnss=xn(I-1, j) - xnu*flux(I,npliq+ks)*dtime
                         if (xnss.lt.0) then
                            dtime1=dtime


                            if (Idiagnose.ge.1) then
                        write(6,'(A,4I5,5E15.6)')'bb solid I-1',
     &      i,ks,j,j1,xn(i-1,j)/xn(I-1,1)*ml(1),dtime,dtime1
                     endif
                     
                        nmin=I
          jmin= j
       isolidmin=ks


            endif
      endif
                 enddo          !K
              endif
              
              enddo
           enddo
c        print*,'dtime d',dtime

c
c     set H2o to maximal 10%
           if (time.lt.timeeq) then
              DO I=1,NS

                 df1 = flux(I+1,1)-flux(I,1)
                 dd=dabs( df1*dtime)/xn(I,1)
c     tolerate 2% change
                 
                 if (dd.ge.0.02) then
                            dtime1=dtime
c     using 10% change to calculate new dtime
                    dtime=0.005/dd * dtime
                    nmin=I
                    jmin=1
                            if (Idiagnose.ge.1) then
                   write(6,'(A,1I5,5E15.6)')'bb H2O I',
     &      i,dd,dtime,dtime1
                   endif
                endif
                 enddo

      if (imode_ph.eq.2) then
c
         DO I =1,NS
            DO J=2, NP-npsolid

                 if (imode_NH4no3.eq.1.and.I.eq.ns) then
                    if (J.eq.18 .or.j.eq.4) then
                    if ( ml(1)*xn(ns,j)/xn(ns,1).le.1D-12) then
                       flux(I+1,j)=flux(I,j)
                    endif
                 endif
                 
                 endif
                 
                 df1 = flux(I+1,j)-flux(I,j)
                 if (df1.ne.0d0 .and. xn(I,j).gt.0d0)then
                 dd=dabs( df1*dtime)/xn(I,j)
c     tolerate 10% change
                 if (dd.ge.0.01) then
                            dtime1=dtime
c     using 10% change to calculate new dtime
                    dtime=0.001/dd * dtime  ! only alow 0.1% change
                    nmin=I
                    jmin=j
                            if (Idiagnose.ge.1) then
                   write(6,'(A,2I5,5E15.6)')'bb I,j',
     &             i,j,dd,dtime,dtime1, df1
c               write(6,'(A,50E16.6)')'flux H2O', (flux(jj,j),jj=2,ns+1)
                   endif
                endif
                endif

               
            enddo
         enddo
      endif
                 
              
              endif
              

c     all major species ( 0.2% of total )  not more than 10%
              
              DO I=1,NS
                     xmc=0d0
                     DO J=6,Npl
                         xmc= xmc+ xn(I,j)
                      enddo

                      DO J=3,npl
                         dfs=0d0
                         DO KS=1,npsolid
                        if (xn(I,ks+npsolid).gt.0d0) then
                    DO K=1,Nsp(KS)
                            j1=nsp_index(KS,K) !species index
                       j0=j1
                       if (j1.eq.25 .or. j1.eq.26 .or.j1.eq.24)j0=24
                       if (j1.eq.29 .or. j1.eq.30 .or.j1.eq.31)j0=29
                       if (iskinetic.eq.0) then
                          if (j1.eq.14 .or.j1.eq.15) J0=5
                       endif
                      if (j0.eq.j)dfs=dfs -xnsp_nv(Ks,K) *flux(I,np+ks)
                      enddo
                   endif
                   
                      
                 if (xn(I+1,ks+npsolid).gt.0d0 .and. I.lt.NS) then
                    DO K=1,Nsp(KS)
                            j1=nsp_index(KS,K) !species index
                       j0=j1
                       if (j1.eq.25 .or. j1.eq.26 .or.j1.eq.24)j0=24
                       if (j1.eq.29 .or. j1.eq.30 .or.j1.eq.31)j0=29
                       if (iskinetic.eq.0) then
                          if (j1.eq.14 .or.j1.eq.15) J0=5
                       endif
                       if (j0.eq.j) 
     & dfs=dfs - xnsp_nv(Ks,K)*flux(I+1,npliq+ks)
 
                         enddo
                      endif
                         

                   enddo        !ks
                   dd=0d0
                   xmj=xn(I,j)*ml(1)/xn(I,1)
                   if (xmj.ge.1D-12)
     &              dd= dabs(flux(I+1,j)-flux(I,j)+dfs)*dtime/xn(I,j)
                   

                    ddmax=.1
c                    if (j.eq.36) dmax=.4
                    ddd=20
                    
c                    if  (xn(I,j)/xmc.le..002d0.and.
c     &           (flux(I+1,j)-flux(I,j)+dfs).gt.0d0) ddmax=.5d0
                    if  (xn(I,j)/xmc.le..002d0.and.
     &           (flux(I+1,j)-flux(I,j)+dfs).gt.0d0) ddmax=.25d0
c     special treatment for calcite dissolution, larger time step 
c                    if (j.eq.14)
c     & print*, 'I, special', I, sshell(I,10),flux(I,NP+10)
c     &            ,       xn(I,46)
                  if (j.eq.14 .and. xn(I,46) .gt. 1D-40 .and.
     % flux(I,NP+10).le.1D-40) then
                       ddmax = .5/ sshell(I,10)
                       
                     ddmax=0.5
                   endif

c     calcite special treament
               if (j.eq.36.and.xn(i,36)/xn(I,1)*ml(1).le.4D-4)ddmax=0.5
              if ( j.eq.36 .and. xn(i,36)/xn(I,1)*ml(1).le.4D-4 .and.
     &                  dfs.lt.0 .and. sshell(i,10).le.10d0
     & .and. (flux(I+1,j)-flux(I,j)+dfs).lt.0d0) then
                      ff= dabs( (flux(I+1,j)-flux(I,j))/dfs)
              if (sshell(i,10).ge.1) ff= (1+(sshell(i,10)-1)/ddd)*ff
                      
c     stationay
c     solid shell i
c              print*,'ff',ff 
                      flux(I,np+10)= ff*flux(I,np+10)
                      flux(I,np+11)= ff*flux(I,np+11)
                      
c     solid shell i+1
                   if (I.lt.NS)flux(I+1,npliq+10)= ff*flux(I+1,npliq+10)
                   if (I.lt.NS)flux(I+1,npliq+11)= ff*flux(I+1,npliq+11)
c     print*,'a'
                      dfs=dfs*ff
                     dd= dabs(flux(I+1,j)-flux(I,j)+dfs)*dtime/xn(I,j)
                    if (Idiagnose.ge.1 ) then
        write(6,'(A,7E15.6)') 'dd', dd, ff, dfs, flux(I+1,j)-flux(I,j)
                     endif
                     ddmax=0.5
                   endif

c     Magnesium
c     calcite special treament
               if (j.eq.32.and.xn(i,32)/xn(I,1)*ml(1).le.4D-4)ddmax=0.5
                   
                   if ( j.eq.32 .and. xn(i,32)/xn(I,1)*ml(1).le.4D-4
     &                   .and.I.le.NS .and.
     &                  dfs.lt.0 .and. sshell(i,12).le.5
     & .and. (flux(I+1,j)-flux(I,j)+dfs).lt.0d0) then
                      ff= dabs( (flux(I+1,j)-flux(I,j))/dfs)

               if (sshell(i,12).ge.1) ff= (1+(sshell(i,12)-1)/ddd)*ff
                      
c     stationay
c     solid shell i
c                      flux(I,np+10)= ff*flux(I,np+10)
                      flux(I,np+12)= ff*flux(I,np+12)
                      
c     solid shell i+1
c                      flux(I+1,np-2)= ff*flux(I+1,np-2)
                   if (I.lt.NS)flux(I+1,npliq+12)= ff*flux(I+1,npliq+12)
c     print*,'a'
                      dfs=dfs*ff
                     dd= dabs(flux(I+1,j)-flux(I,j)+dfs)*dtime/xn(I,j)
                     ddmax=0.5
                   endif

c     do only this is solid flux
                   if (Idiagnose.ge.4) then
                      if (J.eq.24) write(6,'(A,2I5,5E15.6)')
     &                     'bbs', I,J, dd,ddmax,dfs,
     &                      flux(I+1,j)-flux(I,j)+dfs,xn(i,j)
                      
                      
                   endif
                   
                       if ( dd.ge. ddmax .and. dfs.ne.0d0 ) then
                            dtime1=dtime
                            dtime = dtime/dd*ddmax
                    nmin=I
                    jmin=j
                    isolidmin=100
                    if (Idiagnose.ge.1 ) then
      write(6,'(A,2I5,5E15.6)')'bbb species I j',
     &      i,j,dd,dtime,dtime1
      if (j.eq.36) write(6,'(A,I5,16E15.6)')'aa',
     &     I, xn(i,j)/xn(I,1)*ml(1), sshell(I,10),sshell(I,11),
     & flux(I+1,j)-flux(I,j),
     &     dfs ,dd ,ddmax, (flux(I+1,j)-flux(I,j)+dfs)*dtime, xn(i,j) ,
     & (flux(I+1,j)-flux(I,j)+dfs)*dtime+ xn(i,j),
     & (flux(I+1,j)-flux(I,j)+dfs)/(flux(I+1,j)-flux(I,j))

      if (j.eq.36) write(6,'(A,16E15.6)')'bb',
     &   xn(i,NP-2)/xn(I,1)* 55.5,xn(i,NP-1)/xn(I,1)* 55.5
     & ,flux(I, NP+10)      ,flux(I, NP+11) ,flux(I+1, j)-flux(I, j)      
     &, dfs

      if (j.eq.14) write(6,'(A, 18E15.6)') 'm1415',
     &  xn(I,j)/xn(I,1)*ml(1),                   xn(I,15)/xn(I,1)*ml(1)
     & , xn0(I,j)/xn0(I,1)*ml(1),            xn0(I,15)/xn0(I,1)*ml(1)
     & , xn0(I,13)/xn0(I,1)*ml(1),xn(I,j)

      if (j.eq.14) write(6,'(A,6E15.6)')'f co3',
     &     flux(I,np+10), flux(I,np+12),
     &  flux(I+1,j)-flux(I,j),sshell(I,10),sshell(I,12), phshell(I)
      xm6=xn(I,6)*ml(1)/xn(I,1)
      xm7=xn(I,7)*ml(1)/xn(I,1)
      xm13=xn0(I,13)*ml(1)/xn(I,1)
      xm15=xn0(I,15)*ml(1)/xn(I,1)
      xxf=  awshell(I)*xk1t0(I)*xn(I,13)+
     & xk2t0(I)*gamma2(I,7)*xm7* XN(i,13)
      xxb=
     &     xk1tb0(I)*xm6*gamma2(I,6)* gamma2(I,15)* xn0(I,15)
     &     +   xk2tb0(I)*gamma2(I,15)*xn0(I,15)
      if (j.eq.14) write(6,'(A,16E15.6)')'f chm',
     &     xxf,xxb,xxf-xxb,xm6,xm7,xm13,xm15
      if (j.eq.14) write(6,'(A,16E15.6)')'f xkt',
     &     xk1t0(I),xk2t0(I),xk1tb0(I),xk2tb0(I),awshell(I)
                 phh=dlog( xn(I,6)*ml(1)/xn(I,1))/dlog(.1d0)
c       write(6,'(A,16E15.6)')'f gma',
c     &     gamma2(i,6),  gamma2(i,7),  gamma2(i,13),  gamma2(i,15)
c     & ,gamma2(I,14),gamma2(3,6),gamma2(3,14),phh
       if (j.eq.14 .and. gamma2(i,6).ge.1D9)then
         DO jJ=1,NP
            print*,i,jj, xn(I,jj)/xn(I,1)*ml(1), gamma2(i,jj)
            ml(jj)= xn(I,jj)/xn(I,1)*ml(1)
         enddo
         IS=I
         call calHNew(Tdrop,mL)
         call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)

         print*,'gammah > 1D9 , stop',aw
         DO jJ=1,NP
            print*,i,jj, xn(I,jj)/xn(I,1)*ml(1), gamma2(i,jj)
         enddo
         
         DO Ii=1,NS
            print*,ii,jj, xn(Ii,6)/xn(I,1)*ml(1), gamma2(ii,6)
         enddo
         
         stop
      endif
      



                    nmin=i
                    jmin=J
                    
                 endif
                 endif
                         
                      enddo

                      
                  enddo
c     begin integrate
 599              continue
c        print*,'dtime e',dtime

                  xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo
      xnca=0
      DO I=1,ns

         xnca=xnca+ xn(I,36)+xn(I,npl+10)+xn(I,npl+11)         

         enddo

c     end solid
         if (idiagnose.ge.4) print*, 'aa3 begin integr',xnca
         if (idiagnose.ge.4) then

            DO J=1, np
               write(6,'(I5,100E15.6)') J,
     & (xn(I,j), i=1,ns)               
            Enddo
            endif

          if (idiagnose.ge.4) then
               print*,'fluxsolid'
               DO I=1,ns
             DO KS = 1, npsolid
             ISs= NP-NPsolid+KS
             if (flux(i,iss).ne.0d0 .or.flux(i,np+ks).ne.0d0 ) then
           write(6,'(A,2I5,5E15.6)')'aa', i,ks,flux(I,iss),flux(I,np+ks)
             endif
          enddo
          enddo
          
            endif
c        
        
c      flux from solid  to shell i
         DO I=1,ns
          DO KS = 1, npsolid
             ISs= NP-NPsolid+KS

          if (xn(I,NP-npsolid+KS) .gt.0d0) then

           DO K=1,Nsp(KS)
c
                       j1=nsp_index(KS,K) !species index
                       j=j1
                       if (j1.eq.25 .or. j1.eq.26 .or.j1.eq.24)j=24
                       if (j1.eq.29 .or. j1.eq.30 .or.j1.eq.31)j=29
                       if (iskinetic.eq.0) then
                          if (j1.eq.14 .or.j1.eq.15) J=5
                       endif

                       xnu=xnsp_nv(Ks,K) !how many

c     flux from solid  to shell i
            if(j.gt.1) then
                        xn(I, j) =xn(I, j) - xnu*flux(I,np+ks)*dtime
                        endif

                           if (I.ge.2) then
c     flux from solid  to shell i-1
            if(j.gt.1) then
              xn(I-1, j) =     xn(I-1, j) - xnu*flux(I,npliq+ks)*dtime
              endif
              
              endif
                 enddo          !K

c     solidd
                 xn(i,iSs)=xn(I,iSs)+(flux(I,IsS)+flux(I,np+KS))*dtime
c           print*, I, ISS, (flux(I,IsS)+flux(I,np+KS)),xn(i,iss)
           
           if (xn(i,iss).lt.0d0) xn(i,iss)=0d0
          if (xn(i,iss).le.0) then
             print*,'Solid dissolved ', KS
          endif
                 
           endif !xn(I,iss) > 0


        enddo  !KS
        
       enddo !I


          
       
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo
      xnca=0
      DO I=1,ns

         xnca=xnca+ xn(I,36)+xn(I,npl+10)+xn(I,npl+11)         

         enddo
         
         
c     end solid
         if (idiagnose.ge.4) print*, 'aa4 after solid ',dtime
         if (idiagnose.ge.4) then
c            print*,'flux' , flux(5,NP),flux(5,NP+npsolid)

            DO J=1, np
               write(6,'(I5,100E15.6)') J,
     & (xn(I,j), i=1,ns )     , (xn(I,j)*ml(1)/xn(I,1), i=1,ns )      

               
            Enddo
            DO I=1,ns
               DO J=1,np
               ml(j)=ml(1)*xn(I,j)/xn(I,1)
               enddo
            call calHNew(Tdrop,mL)

         call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
         print*,'aw ph', aw, dlog( gammah* ml(6))/dlog(.1d0)
      enddo
      
            
            endif
         

          DO I=1,ns
             vv=4*pi/3 *(x(I+1)**3-x(I)**3)

c             pH= -dlog( ML6shell(I))/dlog(10d0) !pH
             pHa= -dlog( xhshell(I)*ML6shell(I))/dlog(10d0) !pH

             aww=awshell(I)

             if (time.ge.0d0) then
c     mass ratio of liquid
              xm0=xn(I,2)* MM(2)

                DO kk =6,20
                xm0=xm0+ xn(I,kk)* MM(kk)
                enddo


                DO kk =23,np-npsolid
                xm0=xm0+ xn(I,kk)* MM(kk)
                enddo

                
                xm0= xm0/(xm0+ xn(I,1)*mm(1))

                taua=tau_ivea(pha,xm0)

                if (time.le.0) taua=1D10
              xn(I,21)=xn(I,21)  * dexp(-dtime/taua)

             call cal_tau_sars(T,aww,pHa,tauSa)
                if (time.le.0) tausa=1D10
             xn(I,22)=xn(I,22)  * dexp(-dtime/tausa)
          endif
          
          enddo
                
         
 
          DO J=1, Np-npsolid

c            if (j.eq.4) print*, 'j',j, iszero(j)
            if (iszero(j).eq.0 .and. dl_factor2(1,j).gt.0d0) then

                 DO I=1,ns
                  df=(flux(I+1,j)-flux(I,j))
             if (iskinetic.eq.0) then
                if (J.eq.13) df=0d0
                if (J.eq.14) df=0d0
                if (J.eq.15) df=0d0
                if (J.eq.5 .and. I.eq.NS) df=0d0 !equilibrium of CO2 last shell
             endif
             
                if (iskinetic.eq.1) then
                if (J.eq.5  ) df=0d0
                if (J.eq.13 ) df=0d0
                if (J.eq.15 ) df=0d0
                endif
               

              xn(I,j)=xn(I,j)+ df*dtime
              enddo
              
           endif

        enddo


 
c        print*,'dtime f',dtime

       
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo
       

c         end liquid
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo

c     end solid
               xnca=0
      DO I=1,ns

         xnca=xnca+ xn(I,36)+xn(I,npl+10)+xn(I,npl+11)         

         enddo

         if (idiagnose.ge.4) print*, 'aa5 after liq   ',xnca
         if (idiagnose.ge.4) then
            DO J=1, npl 
               write(6,'(I5,100E15.6)') J,
     & (xn(I,j)/xn(I,1)*ml(1), i=1,ns)               
            Enddo
            endif

            
         

            DO I=1,NS
               DO J=1,np-npsolid
            if (xn(I,j).lt.1D-50) then
                 if(j.eq.36 .and. idiagnose.ge.4)
     &                 print*,'liq, Ca', I,xn(I,j)
                 xn(I,j)=0d0
                 
              endif
           enddo
           enddo
           
            if(iskinetic.eq.0) then

c     CO2  equilbrium in shell the outermost shell

       if (partvapco2.ge.1D-40) then
         xkelvin = dexp( 2* sigma * MV(5) /(8.314E7*T*x(NS+1)) )
           pp=partvapco2/xkelvin

c     calculates the CO2, HCO3-, CO2- concentration at given pCO2 partial pressure in the outermost shell

       IS =NS 
       DO J=2,nP
          ML(j)=xn(NS,J)/xn(NS,1)* ml(1)
       enddo

       DO KK=1,NP
       M(KK)=ML(KK)
       enddo
       


       DO KK=1,NP
                ML(kk)=M(KK)
             enddo
             
                HCO2 = 0.034 * dexp(2300 *(1/Tdrop-1/298.15))
                xm13=hco2/1013.5*pp
                ml(13)=xm13
                
                if (Ml(5).le.xm13) ML(5)=xm13



        xk340 = 4.448E-7*dexp(-2133*(1/T-1/298.15)) ! https://www.sciencedirect.
        xk34=xk340/gamma2(NS,6)/gamma2(Ns,15)
          xm15 =xk34*gamma2(NS,1)/(ml(6))

c      HCO3- --> H+ + CO3-2
c     reaction CO2(aq) + H2O --> H+ HCO3-  (1)
c     reaction CO2(aq) + OH- -->   HCO3-   (2)



c     disscociation HCO3- --> H+ + CO3-2

         xk2 = 10D0**(-10.33)*dexp(-3347.3*(1/tdrop-1/298.15d0))
         xk2 = XK2*gamma2(I,15)/gamma2(I,6)/gamma2(I,14)
         xm14= xm15* xk2/ML(6)
        ML(15)= ML(13)*xm15! HCO3-
        ML(14) = ML(13)*xm14    !CO3-2
        ML(5)=ML(13)+ml(14)+ml(15)

        call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)

c        print*,'ppB ', pp, time

       if ( dabs( (pco2-pp)/pp).ge.0.02d0) then
          
          call getm5(Tdrop,ML,pp)
         call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)

         else
        xn(NS,5)= xn(NS,1)*ml(5)/ml(1)
         endif

c        print*,'ppA ', pp, ml(5)

         
       if ( dabs(dlog(pco2/pp)).ge.0.03) then
c          print*, 'before look, after ml5',pco20, pco2
         if(ML(5).le.xm13)ML(5)=xm13
          M(5)=ml(5)
             DO kk=1,np
                ML(kk)=M(kk)
                enddo
                
             ff= pco2/pp
          print*,'ff co2',ff
          
          DO I=1,1000
             if (ff.ge.2) ff=2
             ml(5)= ml(5)*(1 - .001*( ff-1))
       call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
c       print*,i, ml(5), pco2, pco20
       
       if (ff.ge.1 .and. pco2.le. pp) goto 666
       if (ff.lt.1 .and. pco2.ge. pp) goto 666
         
             
             
             
          enddo
          
 666      continue
c          print*, 'after look, after look up'
c                    print*, PCO20*1000,pco2*1000,pp*1000
          

       
      endif
       

        xn(NS,5) =  ML(5)/ML(1)* xn(ns,1)
        xn(NS,13) =  ML(13)/ML(1)* xn(ns,1)
        xn(NS,14) =  ML(14)/ML(1)* xn(ns,1)
        xn(NS,15) =  ML(15)/ML(1)* xn(ns,1)
        
        I=NS
        IS=NS
        
           DO j=1,Np
              ml(j)=ml(1)/xn(I,1)* xn(i,j)
           enddo

           
           call calHNew(Tdrop,mL)

         call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
       phh = dlog(ml(6)*gamma2(I,6))/dlog(.1d0)


        phshell0(i) =phh

        phshell(i) =phh
        awshell(i) =aw
        awshell0(i) =aw

      


      else
        xn(NS,5) =0d0
        xn(NS,13) =0d0
        xn(NS,14) =  0d0
        xn(NS,15) =  0d0
        endif




        endif

             



            if (dtimeb.le.dtime)dtime=dtimeb
            dtimeb=2*dtimeb
c            print*, 'b',time, dtime
c        print*,'dtime bb',dtime
             
cc            if (time+dtime.gt. outputtime(Nout)) 
c     & dtime= outputtime(NoUT)-time+1D-12
            
c        print*,'dtime dd',dtime
 
        
c     NaCL nucleation! keep this part only for EDB experiments
c     for normal runs, it is treated NaCl as all other solid
            
            



c     begin timeloop
            

 
cccccccccccccccc begin nucleation oxalate

 
c      print*,'c finish efflo'
      



       is_Steady=0

cccccccccccccccc
       xn55=0
       xn13s=0
       xsolidc=0
       DO I =1,NS
          xn55=xn55+xn(I,5)
          xn13s=xn13s+xn(I,13)
          xnsolidc=xnsolidc + xn(I,46) + xn(I,48)
       enddo
c

c     begin kinetic
              
         if (idiagnose.ge.6) then
            write(6,'(A,100E15.6)') 'Begin kin. 13', (xn0(i,13), i=1,ns)
            write(6,'(A,100E15.6)') 'Begin kin. 15', (xn0(i,15), i=1,ns)
         endif
                
       imm=0
           if(iskinetic.eq.1) then
 799          continue
            imm=imm+1  
              
              DO I=1,NS
                 xn(I,5)=xn0(I,5)
                 xn(I,13)=xn0(I,13)
                 xn(I,14)=xn0(I,14)
                 xn(I,15)=xn0(I,15)
c                 print*,'ii', xn0(I,15),xn0(i,13)
              enddo
          

       xn13a=0
       xn11=0
       xn15a=0
       DO I=1,NS
          xn13a=xn13a+xn(I,13)
          xn11=xn11+xn(I,1)
          xn15a=xn15a+xn(I,5)
       enddo
        xm13m= xn13a/xn11*ml(1)
        xm15m= xn15a/xn11*ml(1)

        xkelvin = dexp( 2* sigma * MV(1) /(8.314E7*T*x(NS+1)) )
      
      is_steady=1


 
             
          DO I=1,NS


          xk1t =3.7D-2* dexp( -81000./8.314*(1/Tdrop-1/298.15)) !Wang 2010 24.8/10**(-3.7)
          xk1tb =7D4* dexp( -71600./8.314*(1/Tdrop-1/298.15)) !Stumm+Morgan
          xk1tb =12.43D4* dexp( -71600./8.314*(1/Tdrop-1/298.15)) !Wang 2010 24.8/10**(-3.7)
c     reduce xk1tb to match dilted date
c           xk1tb=           xk1tb*.655          
        xk34 = 4.448E-7*dexp(-2133*(1/Tdrop-1/298.15)) ! https://www.sciencedirect.c
        ff= xk34/(xk1t/xk1tb)
          xk1tb =xk1tb/ff !Wang 2010 24.        
c          print*, xk1t/xk1tb, xk34,ff
        
          
          xk2t =12.1E3* dexp( -64000./8.314*(1/Tdrop-1/298.15)) !Wang 2010 24.8/10**(-3.7)
          xk2tb =4E-4* dexp( -114000./8.314*(1/Tdrop-1/298.15)) !Wang 2010 24.8/10**(-3.7)

          xkw=  dexp( -0.92644d1-0.68727E+04/tdrop) !dissociation of H2O
          xk11= xkw* xk2t/xk2tb
          ff= xk34/xk11
          xk2tb =4E-4* dexp( -114000./8.314*(1/Tdrop-1/298.15))/ff !Wang 2010 24.8/          
          xk11= xkw* xk2t/xk2tb


           xm15= xn(i,14)/xn(I,1)*ml(1) + xn(i,26)/xn(I,1)*ml(1)
           xm15= xm15+xn(i,28)/xn(I,1)*ml(1)
           xm15= xm15+xn(i,31)/xn(I,1)*ml(1)

           
           xmcl=0
           DO J=8,np
              if (izc(j).ge.-1) xmcl=xmcl+ xn(I,j)/xn(I,1)*ml(1)
           enddo
           
          xm6=xn(i,6)/xn(I,1)*ml(1)
          xm7=xn(I,7)/xn(I,1)* ml(1)
          phh= dlog(xm6*gamma2(I,6))/dlog(.1d0)

          xcatkm = 10d0**( catkm +dcatkm*(phh-7d0))

          ml(35)= xn(I,35)*ml(1)/xn(I,1)
           fenzyme= xcatkm/xk1t*ml(35)
           ff= 1+fenzyme+fCO3*xm15+fcl*xmcl

           xk2tb0(I)=xk2tb*ff
           xk2t0(I)=xk2t*gamma2(I,13)*ff
           xk1tb0(I)=xk1tb*ff
           xk1t0(I)=xk1t*ff*gamma2(I,13)

          freact = awshell(I)*xk1t*xn(I,13)+ !CO2 + H2O -- H+ + HCO3-
     & xk2t *  xm7 * xn(1,13)*gamma2(I,7)  ! CO2 + OH1
     &         - xn(I,15) *
     &         (gamma2(I,6)*gamma2(I,15)*xm6*xk1tb+xk2tb*gamma2(I,15) ) ! HCO3- + H+ and HCO3- -> OH- +

          
c     3      - + H+ and HCO3- -> OH- +CO2
          
        enddo

cccccccccccccccc
       if (partvapco2.ge.1D-40 .or. xn(1,5).ge.1D-50) then
          
          DO I=1,ns
             if (xn(I,5).le.0d0) then
                print*, 'CO2 molality is zero with CO2 in the gas phase'
                print*, 'If you really want to study CO2 uptake, '
        print*, 'input a minor concentration of CO2 5. line in ML.dat'
             stop
             endif
             
          enddo
          

          xkelvin = dexp( 2* sigma * MV(5) /(8.314E7*T*x(NS+1)) )
           pp=partvapco2/xkelvin

c     calculates the CO2, HCO3-, CO2- concentration at given pCO2 partial pressure in the outermost shell

       DO J=2,nP
          ML(j)=xn(NS,J)/xn(NS,1)* ml(1)
       enddo

       IS =NS 
       
          xkelvin = dexp( 2* sigma * MV(5) /(8.314E7*T*x(NS+1)) )
           pp=partvapco2/xkelvin

           call calhnew(Tdrop,ml)
           
          HCO2 = 0.034 * dexp(2300 *(1/Tdrop-1/298.15))/gamma2(ns,13)
           
       
c
c     calculates the CO2, HCO3-, CO2- concentration at given pCO2 partial pressure in the outermost shell
c      if (NS.eq.I) then
c     call aw_back(tdrop,ML,aw0,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      
      icc=0
      

c      if (is_steady.eq.1) then
         if (NS.gt.50) then
          print*, 'error: increase dimension of matrix!'
          stop

       endif
 
c
         DO I=1,ns
         vshellliq(I )= xn(I,1)*mv(1)+xn(I,2)*mv(2)
         DO j=6,Np-NPsolid
            vshellliq(I )= vshellliq(I )+xn(I,j)*mv(j)
            
         enddo
      enddo
      
      
      DO I=1,ns
          xn130(I)=xn(I,13)
          xn13(I)=xn(I,13)
          a15(I)= xn(I,15)*ml(1)/xn(I,1) * gamma2(I,15)
          c15(I)= xn(I,15)/vshellliq(I)
          a13(I)= xn(I,13)*ml(1)/xn(I,1) * gamma2(I,13)
          c13(I)= xn(I,13)/vshellliq(I)

       enddo

       DO I=1+NS,ns+NS
          II=I-NS
          xn130(I)=xn(Ii,15)
          xn13(I)=xn(iI,15)
c          print*, 'I', I,xn(Ii,15)

       enddo
       
c       xk1t0=xk1t
c       xk2t0=xk2t
       
       
 245   continue
       if (NS.ge.2) then
       icc=icc+1
      J=13
      DH2O0=0.211*1013.d0/PRESS*(T/273.15d0)**1.94 ! Pruppacher + Klett
      velocity=(8.d0*8.314D7*T/PI/MM(1))**0.5
      rad=x(NS+1)
      if (venti.le.1) venti=1
      
      DH2O=venti*DH2O0/(1d0 +
     +        (4.d0*DH2O0/RAD/velocity))
         DGas =DH2O* dsqrt(MM(1)/MM(J))
         AA= 4*pi * x(NS+1) *Dgas/8.314/T * 1013.5E-4

         
         DO I=1,ns
            
                  if (i.le.NS-1) then
               aw= (awshell(I)+awshell(I+1))/2
c     abcd
c               xvv=mv(j)
               DO k=1,NP

              ml0(k)=ml(1)/2d0*(xn(I+1,k)/xn(I+1,1)+xn(I,k)/xn(I,1))

            enddo
         call caldl(tdrop,aw,ml0,x,J,dl)
            
         dl=dl_factor2(1,J)*dl ! diffusion coefficient
         dln(I)=  dl  !interface I,I+1
         deltaxn(I)= (x(I+2)-x(I))/2 
            endif
            
         enddo

         J=15
         DO I=1+NS,ns+NS
         II=I-NS
         vshellliq(I )= vshellliq(Ii)
         if (i.le.2*NS-1) then
               aw= (awshell(Ii)+awshell(Ii+1))/2
            ml0(k)=ml(1)/2d0*(xn(Ii+1,k)/xn(Ii+1,1)+xn(iI,k)/xn(iI,1))
            call caldl(tdrop,aw,ml0,x,J,dl)
           
               dl=dl_factor2(1,15)*dl ! diffusion coefficient
         dln(I)=  dl  !interface I,I+1
         deltaxn(I)= (x(Ii+2)-x(Ii))/2 
            endif
            
         enddo

         
             DO I=1,NS
c     calculte matrix

            xk1t=xk1t0(I)
            xk2t=xk2t0(I)
            xk1tb=xk1tb0(I)
            xk2tb=xk2tb0(I)
            
            DO j=1,2*NS
                   matrix(I,j)=0
                enddo
c     define matrix liquid phase diffusion
                
                
c     chemisty
                xm7= xn(I,7)*ml(1)/xn(I,1)
                xm6= xn(I,6)*ml(1)/xn(I,1)

            matrixf(I)= xn13(I) - xn130(i)

c     flux from i to i-1
            if(I.ge.2) then
c        gx=(xn13(i-1)/vshellliq(i-1)-xn13(I)/vshellliq(i))/deltaxn(I-1)
        gx=cm2(I,13)* dlog(A13(I-1)/A13(I))/deltaxn(I-1)
        
        matrixf(I)=matrixf(I)      -dtime*4*pi*
     &gx*x(I)**2*dln(i-1)
!from i-1 to i
            endif

            if(I.lt.NS) then
c        gx=(xn13(i+1)/vshellliq(i+1)-xn13(I)/vshellliq(i))/deltaxn(I)
        gx=cm2(I+1,13)* dlog(A13(I+1)/A13(I))/deltaxn(I)
         matrixf(I)=matrixf(I)      -dtime*4*pi*
     &x(I+1)**2* gx*dln(i) !from i+1 to i
      endif

c     gas phase

      
      if(I.eq.ns)then
               matrixF(NS)=matrixF(ns)+
     & DTIME*aA*mL(1)/HCO2/xn(I,1)*xkelvin*xn13(I)
     & -DTIME*AA* partvapco2/1013.5d0 
      endif

C     CHEMISTRY
            matrixf(i)=matrixF(i)-DTIME*
     &      (-awshell(I)*xk1t*xn13(I)-xk2t*gamma2(I,7)*xm7* XN13(i)
     & +xk1tb*xm6*gamma2(I,6)* gamma2(I,15)* xn13(I+NS)
     &     +   xk2tb*gamma2(I,15)*xn13(I+NS))
C     FINITO  f
            
c            define derivative matrix
c     chemistRY, LIQUID, GAS
c     dyi/dxi
c     f(I,j)=1 - dtime * dfi/dxi =0 
            if (I.lt.NS .and. I.ge.2) then
c               ff1=-1/vshellliq(I)
c               ff2=-1/vshellliq(I)
               ff1= 1/2d0*dlog(a13(I-1)/A13(I))/vshellliq(I)
     &              - cm2(I,13) /xn13(I)
               ff2= 1/2d0*dlog(a13(I+1)/A13(I))/vshellliq(I)
     &              - cm2(I+1,13) /xn13(I)
               
           matrix(I,I) = 1-dtime*(-awshell(I)*xk1t-xk2t*gamma2(I,7)*xm7) !chem
     & -DTIME*(4*PI*x(I)**2*ff1/deltaxn(I-1)*dln(i-1)) !from i to i-1
     & -DTIME*(4*PI*x(I+1)**2*ff2/deltaxn(I)*dln(i)) !from i to i+1
        endif
        
            if ( I.eq.1) then
               ff2=-1/vshellliq(I)
               ff2= 1/2d0*dlog(a13(I+1)/A13(I))/vshellliq(I)
     &              - cm2(I+1,13) /xn13(I)

               
          matrix(I,I) = 1-dtime*(-awshell(I)*xk1t-xk2t*gamma2(I,7)*xm7) !chem
     & -DTIME*(4*PI*x(I+1)**2*ff2/deltaxn(I)*dln(i)) !from i to i+1
      endif
      if (I.eq.NS ) then
               ff1=-1/vshellliq(I)
               ff1= 1/2d0*dlog(a13(I-1)/A13(I))/vshellliq(I)
     &              - cm2(I,13) /xn13(I)
         
       matrix(I,I) = 1-dtime*(-awshell(I)*xk1t-xk2t*gamma2(I,7)*xm7) !chem
     & -DTIME*(4*PI*x(I)**2*ff1/deltaxn(I-1)*dln(i-1)) !from i to i-1
     & +DTIME*aA*mL(1)/HCO2/xn(I,1)*xkelvin
      endif

      
      
      
            
c     liquid phase diffusion
c     f(I,j)= - dtime * dfi/dxj 

                if (I.ge.2) then
c     flux I-1 -_> I
                   ff1=1/vshellliq(I-1)
               ff1= 1/2d0*dlog(a13(I-1)/A13(I))/vshellliq(I-1)
     &              + cm2(I,13) /xn13(I-1)
                   matrix(I,I-1) = 
     & -DTIME*( 4*PI*x(I)**2*ff1/deltaxn(I-1)*dln(i-1)) !from i-1 to i
      endif

                if (I.lt.NS) then
                   ff2=1/vshellliq(I+1)
               ff2= 1/2d0*dlog(a13(I+1)/A13(I))/vshellliq(I+1)
     &              + cm2(I+1,13) /xn13(I+1)

                   matrix(I,I+1) = 
     & -DTIME*( 4*PI*x(I+1)**2*ff2/deltaxn(i)*dln(i)) !from i-1 to
                endif

c     chemsity HCO3
                matrix(i,I+NS)=-DTIME*(
     & +xk1tb*xm6*gamma2(I,6)* gamma2(I,15)
     &     +   xk2tb*gamma2(I,15))
                
               enddo

               
c     HCO3-
               DO I=1+NS,NS+NS
                  II=I-NS
c     calculte matrix
            xk1t=xk1t0(Ii)
            xk2t=xk2t0(Ii)
            xk1tb=xk1tb0(Ii)
            xk2tb=xk2tb0(Ii)


                  DO j=1,2*NS
                   matrix(I,j)=0
                enddo
c     define matrix liquid phase diffusion
                
                
c     chemisty
                xm7= xn(Ii,7)*ml(1)/xn(Ii,1)
                xm6= xn(Ii,6)*ml(1)/xn(Ii,1)

            matrixf(I)= xn13(I) - xn130(i)

c     flux from i to i-1
            direction=1d0
            if(II.ge.2) then

c         gx=(xn13(i-1)/vshellliq(i-1)-xn13(I)/vshellliq(i))/deltaxn(I-1)
        gx=cm2(II,15)* dlog(A15(II-1)/A15(II))/deltaxn(I-1)
c        print*,'gg', gx ,gx1 
        matrixf(I)=matrixf(I)      -dtime*4*pi*
     &gx*x(Ii)**2*dln(i-1)
       matrixf(I)=matrixf(I) -dtime * cm2(II,15)*velarr(Ii-1)
     &       * izc(15)*(-1d0)*direction        

!     from i-1 to i
            endif

            if(II.lt.NS) then
c        gx=(xn13(i+1)/vshellliq(i+1)-xn13(I)/vshellliq(i))/deltaxn(I)
        gx=cm2(II+1,15)* dlog(A15(II+1)/A15(II))/deltaxn(I)
c        print*,'gg', gx ,gx1 

        matrixf(I)=matrixf(I)      -dtime*4*pi*
     &x(Ii+1)**2* gx*dln(i) !from i+1 to i
        matrixf(I)=matrixf(I) -dtime * cm2(II+1,15)*velarr(Ii)
     &       * izc(15)*direction        

      endif


C     CHEMISTRY
            matrixf(i)=matrixF(i)-DTIME*
     &      (awshell(Ii)*xk1t*xn13(iI)+xk2t*gamma2(iI,7)*xm7* XN13(ii)
     & -xk1tb*xm6*gamma2(Ii,6)* gamma2(Ii,15)* xn13(I)
     &     -   xk2tb*gamma2(iI,15)*xn13(I))
C     FINITO  f
            
c            define derivative matrix
c     chemistRY, LIQUID, GAS
c     dyi/dxi
c     f(I,j)=1 - dtime * dfi/dxi =0 
            if (iI.lt.NS .and. iI.ge.2) then
c               ff1= -1/vshellliq(I)
c               print*,ii, a15(II-1), a15(II),vshellliq(I)
               
               ff1= 1/2d0*dlog(a15(II-1)/A15(II))/vshellliq(I)
     &              - cm2(II,15) /xn13(I)
               
               
               ff2= -1/vshellliq(I)
               ff2= 1/2d0*dlog(a15(II+1)/A15(II))/vshellliq(I)
     &              - cm2(II+1,15) /xn13(I)

c     &-1/2d0/vshellliq(I)-
c     &              xn13(I+1)/xn13(I)/vshellliq(I+1)/2



               
               matrix(I,I) = 1-dtime*( !chem
     & -xk1tb*xm6*gamma2(Ii,6)* gamma2(Ii,15)
     &     -   xk2tb*gamma2(iI,15))
     & -DTIME*(4*PI*x(Ii)**2*ff1/deltaxn(I-1)*dln(i-1)) !from i to i-1
     & -DTIME*(4*PI*x(Ii+1)**2*ff2/deltaxn(I)*dln(i)) !from i to i+1

c     vel  II.eq.NS
               matrix(I,I)
     &              =matrix(I,I)-dtime*1/2d0/vshellliq(I)*izc(15)
     & * direction*velarr(iI-1)*(-1d0)
        matrix(I,I)= matrix(I,I)-dtime*1/2d0/vshellliq(I)*izc(15)
     & * direction*velarr(Ii)

            endif

            if ( iI.eq.1) then
                               ff2= -1/vshellliq(I)
               ff2= 1/2d0*dlog(a15(II+1)/A15(II))/vshellliq(I)
     &              - cm2(II+1,15) /xn13(I)

c     &-1/2d0/vshellliq(I)-
c     &              xn13(I+1)/xn13(I)/vshellliq(I+1)/2

                   

                           matrix(I,I) = 1-dtime*( !chem
     & -xk1tb*xm6*gamma2(Ii,6)* gamma2(Ii,15)
     &     -   xk2tb*gamma2(iI,15))
     & -DTIME*(4*PI*x(Ii+1)**2*ff2/deltaxn(I)*dln(i)) !from i to i+1

c     vel  II.eq.NS
               matrix(I,I)
     &              =matrix(I,I)-dtime*1/2d0/vshellliq(I)*izc(15)
     & * direction*velarr(iI)

                        endif
            if (Ii.eq.NS ) then
                              ff1= -1/vshellliq(I)
               ff1 = 1/2d0*dlog(a15(II-1)/A15(II))/vshellliq(I)
     &             -  cm2(II,15) /xn13(I)

c     &-1/2d0/vshellliq(I)-
c     &              xn13(I-1)/xn13(I)/vshellliq(I-1)/2



               matrix(I,I) = 1-dtime*( !chem
     & -xk1tb*xm6*gamma2(Ii,6)* gamma2(Ii,15)
     &     -   xk2tb*gamma2(iI,15))
     & -DTIME*(4*PI*x(iI)**2*FF1/deltaxn(I-1)*dln(i-1)) !from i to i-1
c     vel  II.eq.NS
               matrix(I,I)
     &              =matrix(I,I)-dtime*1/2d0/vshellliq(I)*izc(15)*(-1)
     & * direction*velarr(iI-1)

c               PRint*, 'NS', i, vshellliq(I),deltaxn(I-1)
c               print*, 'ns', matrix(i,i)

            endif

            

c     f(I,j)= - dtime * dfi/dxj 

                if (iI.ge.2) then
                   ff1=1/vshellliq(I-1)
                   ff1= 1/2d0*dlog(a15(II-1)/A15(II))/vshellliq(I-1)
     &              + cm2(II,15) /xn13(I-1)

c     &+1/2d0/vshellliq(I-1)+
c     &              xn13(I)/xn13(I-1)/vshellliq(I)/2

                   
            matrix(I,I-1) = 
     & -DTIME*( 4*PI*x(Ii)**2*ff1/deltaxn(I-1)*dln(i-1)) !from i-1 to i
c     vel
               matrix(I,I-1)
     &            =matrix(I,I-1)-dtime*1/2d0/vshellliq(I-1)*izc(15)*(-1)
     & * direction*velarr(iI-1)


         endif

                if (Ii.lt.NS) then
                   ff2=1/vshellliq(I+1)
                   ff2= 1/2d0*dlog(a15(II+1)/A15(II))/vshellliq(I+1)
     &              + cm2(II+1,15) /xn13(I+1)

c     &                  +1/2d0/vshellliq(I+1)+
c     &              xn13(I)/xn13(I+1)/vshellliq(I)/2


                   matrix(I,I+1) = 
     & -DTIME*( 4*PI*x(iI+1)**2*ff2/deltaxn(i)*dln(i)) !from i-1 to

            
c     vel
               matrix(I,I+1)
     &           =matrix(I,I+1)-dtime*1/2d0/vshellliq(I+1)*izc(15)
     & * direction*velarr(iI)


                endif

c     chemsity HCO3
                matrix(i,Ii)=-DTIME*
     &      (awshell(Ii)*xk1t +xk2t*gamma2(iI,7)*xm7)

               enddo

            endif               !nS >=2
c     NS.eq.1
            if (NS.eq.1) then


       icc=icc+1
      J=13
      DH2O0=0.211*1013.d0/PRESS*(T/273.15d0)**1.94 ! Pruppacher + Klett
      velocity=(8.d0*8.314D7*T/PI/MM(1))**0.5
      rad=x(NS+1)
      DH2O=DH2O0/(1d0 +
     +        (4.d0*DH2O0/RAD/velocity))
         DGas =DH2O* dsqrt(MM(1)/MM(J))
         AA= 4*pi * x(NS+1) *Dgas/8.314/T * 1013.5E-4

         
         
             DO I=1,NS
c     calculte matrix
            xk1t=xk1t0(I)
            xk2t=xk2t0(I)
            xk1tb=xk1tb0(I)
            xk2tb=xk2tb0(I)
c            xk1t=xk1t0*gamma2(II,13)
c            xk2t=xk2t0*gamma2(II,13)
           

                DO j=1,2*NS
                   matrix(I,j)=0
                   matrix(I+ns,j)=0
                enddo
c     define matrix liquid pha
                
                
c     chemisty
                xm7= xn(I,7)*ml(1)/xn(I,1)
                xm6= xn(I,6)*ml(1)/xn(I,1)

            matrixf(I)= xn13(I) - xn130(i)

c     gas phase

               matrixF(NS)=matrixF(ns)+
     & DTIME*aA*mL(1)/HCO2/xn(I,1)*xkelvin*xn13(I)
     & -DTIME*AA* partvapco2/1013.5d0 


C     CHEMISTRY
            matrixf(i)=matrixF(i)-DTIME*
     &      (-awshell(I)*xk1t*xn13(I)-xk2t*gamma2(I,7)*xm7* XN13(i)
     & +xk1tb*xm6*gamma2(I,6)* gamma2(I,15)* xn13(I+NS)
     &     +   xk2tb*gamma2(I,15)*xn13(I+NS))
            

       matrix(I,I) = 1-dtime*(-awshell(I)*xk1t-xk2t*gamma2(I,7)*xm7) !chem
     & +DTIME*aA*mL(1)/HCO2/xn(I,1)*xkelvin


          
c     chemsity HCO3
                matrix(i,I+NS)=-DTIME*(
     & +xk1tb*xm6*gamma2(I,6)* gamma2(I,15)
     &     +   xk2tb*gamma2(I,15))
                
               enddo

               
c     HCO3-
               DO I=1+NS,NS+NS
                  II=I-NS
c            xk1t=xk1t0*gamma2(iI,13)
c            xk2t=xk2t0*gamma2(Ii,13)
            xk1t=xk1t0(Ii)
            xk2t=xk2t0(Ii)
            xk1tb=xk1tb0(Ii)
            xk2tb=xk2tb0(Ii)

                  
c     calculte matrix
                DO j=1,2*NS
                   matrix(I,j)=0
                enddo

                
                
c     chemistry
                xm7= xn(Ii,7)*ml(1)/xn(Ii,1)
                xm6= xn(Ii,6)*ml(1)/xn(Ii,1)

            matrixf(I)= xn13(I) - xn130(i)

c     flux from i to i-1


C     CHEMISTRY
            matrixf(i)=matrixF(i)-DTIME*
     &      (awshell(Ii)*xk1t*xn13(iI)+xk2t*gamma2(iI,7)*xm7* XN13(ii)
     & -xk1tb*xm6*gamma2(Ii,6)* gamma2(Ii,15)* xn13(I)
     &     -   xk2tb*gamma2(iI,15)*xn13(I))
               
               matrix(I,I) = 1-dtime*( !chem
     & -xk1tb*xm6*gamma2(Ii,6)* gamma2(Ii,15)
     &     -   xk2tb*gamma2(iI,15))
      
            
                matrix(i,Ii)=-DTIME*
     &      (awshell(Ii)*xk1t +xk2t*gamma2(iI,7)*xm7)

               enddo

            endif               !nS =1

            
            
                  ns2=2*ns
c                  print*,'bef ma', matrix(1,1)
                call matrix_inv(matrix,matrixinv,matrixb,x13s,ns2)
c                  print*,'af ma', matrixinv(1,1),matrixinv(ns2,ns2)
       
        irepeatcond=0
        dd=0d0

       DO I=1,2*ns
      xx=0
      DO J=1,2*ns
         xx=xx+ matrixinv(I,J)* matrixf(j)
      enddo
      
     
      xn13(I)= xn13(I) - xx


      

             if (dabs(xn13(I))/xn(I,1)*ml(1).ge. 1D-30 ) then

                if ( dabs( xx/xn13(I)).ge.1D-4) then
                 irepeatcond = 1
                 dd= dabs( xx/xn13(I))

                 II=I
             if(I.gt.nS) II=I-ns
       if (icc.ge.10 .and. idiagnose.ge.1) write(6,'(A,2I5,60E15.6)')
     &                 'matrix, converg', I,icc, xx/xn(Ii,1)*ml(1),
     &                 xn13(I)/xn(Ii,1)*ml(1),dd,dtime

       
              endif

              if (icc.ge.10) then
                 dtime=dtime/5
                 goto 799
              endif
              
                endif

      enddo

c     calculate the new gammas
         isczero=0
      
          DO I =1,NS
c               print*,i, xm13
c             print*, 'I13', i,xn13(I), xn(I,13)
c             print*, 'I15', i, xn13(I+NS), xn(I,15)
             
           xn(I,13)=xn13(I)
           xn(I,15)=xn13(I+NS)
           xn(I,5) = xn(I,13)+xn(I,14)+ xn(i,15)
           if (xn13(I).le. 0d0) isczero=1
           if (xn13(I+ns).le. 0d0) isczero=1
           
           DO J=1,NP-npsolid
              ml(J)= xn(I,j)/xn(I,1)*ml(1)
          enddo
           IS=I
c     no new equilibrium calculation and gamma updates
           
c     call calhnew(Tdrop, ml)
C     IGNORE THE CHANGE OF VOLUME
           a15(i) = XN(i,15)/XN(i,1) *  ML(1)*GAMMA2(i,15)
           c15(I)= xn(I,15)/vshellliq(I)
        enddo
        if ( isczero.eq.1) then
           dtime=dtime/5
          if (idiagnose.ge.1) print*,' n13,15 become 0, reduce dtime'
     & ,imm,dtime
           goto 799
        endif
        
        
      if (       irepeatcond.eq.1 .and. icc.le.50) goto 245

      if (       irepeatcond.eq.0 ) goto 247

        if (Icc.ge.50) then
           print*, 'warning: Newton raphson not  converge ',icc,
     &           irepeatcond,dd
         endif
 247     continue
         
       DO I =1,NS
c               print*,i, xm13

           xn(I,13)=xn13(I)
           xn(I,15)=xn13(I+NS)
           xn(I,5) = xn(I,13)+xn(I,14)+ xn(i,15)
               
            enddo
       
c            print*, 'dtime, pH', time, dtime,
c     & dlog(ml(6)*gamma2(1,6))/dlog(.1d0)
            
               
             xm13= xn(NS,13)*ML(1)/xn(NS,1)
         

             
         endif  !end pco2=0

        
      
      I=NS

ccccccccccccccc
        


        I=NS
      
      Do J =1, NP
         ml(J)= xn(I,j)/xN(I,1)*ml(1)
      enddo
      is=NS
      
      pco2q = ml(13) * 1013.5 /HCO2
      call vapnew(tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
      
      endif !end iskinetic.eq.1
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

         enddo
         xnca=0d0
         
      DO I=1,ns

         xnca=xnca+ xn(I,36)+xn(I,npl+10)+xn(I,npl+11)         

         enddo

         if (idiagnose.ge.4) print*,'aa6, after kinet',xnca
         if (idiagnose.ge.4) then
            DO J=1, npl 
               write(6,'(I5,100E15.6)') J,
     & (xn(I,j)/xn(I,1)*ml(1), i=1,ns)               
            Enddo

         endif

          
         if (idiagnose.ge.6) then
            write(6,'(A,100E15.6)') 'After kin. 13', (xn(i,13), i=1,ns)
            write(6,'(A,100E15.6)') 'After kin. 15', (xn(i,15), i=1,ns)
         endif


         if (iskinetic.eq.1 ) then
          DO I=1,NS
             xn(I,5)= xn(I,13)+xn(I,14)+xn(I,15)
             xn(i,14)= xn(I,14)+xn(I,15)
             xn(I,15)=0d0
          enddo
         endif
                
       
c     finish CO2 kinetic
c       print*, xn55, xn13s
       
c       if (xn55.le. xn13s*2 .and. xnsolidc.gt.0d0) iskinetic=0
c       if (xn55.le. xn13s*2 .and. xnsolidc.gt.0d0) iskinetic0=0
       
           Tdrop=ta              ! the droplet T is calculated in cal_flux

           if (isdiv.eq.2) goto 699
           
         if (partvapco2.ge.1D-40 .or. xn(NS,5).ge.1D-50) then
         if (iskinetic.ge.1) then
       DO I=1,NS
          if (xn(I,14)+xn(I,15).le.0d0) then
             
             DO j=1,NP
                print*,'j', j, xn(I,j)*ml(1)/xn(I,1)
             enddo
             print*,'error: carbonate ions <= 0'
             stop
          endif
       enddo
      endif
      endif
      
       
c      print*,'dtime cc',dtime
      

      pco2eq=pco2               ! for output
      dphmax=0d0
      dawmax=0d0

      do I=1,ns
         DO J=1,npliq
            m(J)= m(1)*xn0(I,j)/xn0(I,1)
            ml(J)= ml(1)*xn(I,j)/xn(I,1)
            gamma2(i,j)= gamma20(i,j)
         enddo
         IS=I

      if ( iseqAA.eq.1.and.I.eq.NS) then
         m(3)=ml(3)
      endif
      if ( iseqNH3.eq.1 .and. I.eq.NS) then
         m(4)=ml(4)
      endif
         
       DO J=1,Npl
c                    gamma2(i,j)= gamma20(i,j)
       enddo
         call calHNew(Tdrop,m)
          call aw_back
     & (tdrop,M,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
       phshell0(I) = dlog(m(6)*gamma2(I,6))/dlog(.1d0)
       awshell0(I)= aw
       phshell(I) = phshell0(i)
       awshell(I)= aw
       gammah0=gammah
         IS=I

         
        call calHNew(Tdrop,mL)
          call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
       phh = dlog(ml(6)*gamma2(I,6))/dlog(.1d0)
c      if (idiagnose.ge.1.and.i.eq.nmin) print*, 'phh', phshell0(i), phh
       dphh= dabs( phh-phshell0(I))
       daw= dabs( aw-awshell0(I))
c     set back the activity coefficients
       
       if (kka.ge.99) then
          DO J=1,npliq
             if (ml(j).gt.0d0 .or. m(j).gt.0d0)
     &        write(6,'(A, 2I5,5E15.6)') 'kka ',i,j, ml(j), ml(j)-m(j)

          enddo
        write(6,'(A,5E15.6)') 'ph0,ph ',phshell0(i), phh,gammah0,gammah

       endif
c     no pH check when pHmode=2
       
       if (imode_ph.eq.2) dphh=0
       
       if (iskinetic.eq.0d0 .and. I.eq.NS) dphh=0d0
       if (dphh.gt. dphmax) then
        dphmax= dphh
        phhmax=phh
        imax=I
c        jmin=2001
      endif

      if (daw.gt. dawmax .and. (time.le.timeeq) ) then
         dawmax= daw
         awmax=aw
        imaxaw=I
c        nmin=I
c        jmin=1001
      endif


      enddo


      if (idiagnose.ge.1) then

c         write(6,'(A,3I5,14E15.6)') 'dpH',kka,Imax,isdiv,
c     & dphmax,phhmax,phshell0(Imax)


      
      if (gamma2(imax,6).ge.1D6) then
      DO jj=1,np
         ml(jj)= xn(imax,jj) *ml(1)/xn(imax,1)
         print*, jj, ml(jj), gamma2(imax,jj)
      enddo
      IS=Imax
      print*,'dddd'
      call calhnew(Tdrop, ml)
          call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
       phh = dlog(ml(6)*gamma2(Imax,6))/dlog(.1d0)
      DO jj=1,np

         print*, jj, ml(jj), gamma2(imax,jj)
      enddo

      print*, 'aw', aw,phh
      
       stop
      
      endif
      endif
      
      

      if (dphmax.ge. 0.05d0) then

         if (kka.ge.10) then
           write(6,'(A,5E15.6)') 'kka', dphmax,dawmax,dtime
            endif
         Jmin=1001


         isdiv=1
         if (phhmax.le.10)then
      dtime=dtime/10d0
      nmin=imax
      jmin=1001
      
      goto 499

      endif
      
         dtime=dtime/dphmax/2*.05d0
      nmin=imax
      jmin=1001
         goto 499
      endif
      dawmax0=0.01
      if (i.ge.ishelleq)dawmax0=0.05
      
      if (dawmax.ge. dawmax0 .and.I.lt.ishelleq) then
        Jmin=101
      nmin=imaxaw
      dtime1=dtime
         
         dtime=dtime/dawmax/2*dawmax0
         if (Idiagnose.ge.1) write(6,'(A,I5,5E15.6)') 'aw dtime', Nmin,
     & dtime1,dtime         
         isdiv=1
         goto 499
      endif
c      print*,'dtime ee',dtime
      
      
      if (dtime.gt.10.0001d0 ) then
         dtime=10d0
            isdiv=2
            goto 499
      endif

      if (dtimeb.le.dtime)then
           dtime=dtimeb
            isdiv=2
           
            dtimeb=2*dtimeb
            goto 499
            endif
            
c      print*,'dtime ff',dtime

c      if (time+dtime.gt. outputtime(Nout)) then
c       dtime= outputtime(NoUT)-time+1D-12
       
c       isdiv=2
c            goto 499
c                       endif
                       
      
c     else
 699                   continue
                       

      if(isdiv.eq.0 ) then
         if(dphmax.le..01 .and. time.ge.timeeq) then
            dtime=dtime*2
            if (Idiagnose.ge.1) then
            print*, ' increase dtime timeeq', dtime
            endif
            isdiv=1
            
            goto 499
         endif
c     increase dtime
         if(dphmax.le..025 .and. dawmax.le.5D-3) then
            dtime=dtime*2
            if (Idiagnose.ge.1) then
            print*, ' increase dtime', dtime
            endif
            isdiv=1
            
            goto 499

         endif


c     increase dtime
      ENDIF


      
 345  continue

       dtime00=dtime                
      
       time=time+dtime
       times=times+dtime
       dtimesplit=dtimesplit+dtime
       sdtime=sdtime+dtime

      
c
c      print*,dtimesplit, dtimesplitmax,imode_NH4NO3
c     Set aw equilibrium

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print*,'IshellEQ', ishelleq
      
       DO I=Ishelleq,NS
            aw=awshell(I)

         isH2O=1

            if (ish2o.eq.1) then
c     check aw
             xkelvin = dexp( 2* sigma * 18d0 /(8.314E7*T*x(NS+1)) )
         rhh= rh/xkelvin
c c        print*,'rhh',rhh,idiagnose
         
         
c         dd=.3D-3
          dd=2D-3
c         print*,aw,rhh
c          print*,'aww', i, aw, aw- rhh

        if (dabs(Aw-rhh).ge.dd) then

           IS=I
           aw1=aw*xkelvin
           DO j=2,Np
              ml(j)=ml(1)*xn(I,j)/xn(I,1)
           enddo
           ml0=ml
           call cal_ml(Ta,rhh,ML)
           xnI1=xn(I,1)
        xn(I,1) = (ML0(2)+ml0(16)+ML0(3)+ml0(5)+ml0(18)+ml0(21))/
     &          (ML(2)+ml(16)+ML(3)+ml(5)+ml(18)+ml(21)) *xni1

c           xn(I,1) = (ML0(2))/
c     &          (ML(2)) *xni1
           

           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           awshell(I)=Aw
c        print*,'xn1', xn(I,1),xni1,i
        
c      repeat check
          ff= 1000d0 /(MM(1)*xn(I,1))

        
          DO jJ=2,NP
             xx=ml(jj)
             ML(jJ)= xn(I,jj)*ff
             
c             if (ml0(jj).gt.0d0) write(6,'(I5,5E15.6)')
c     &             jJ,ml0(jj),xx/ml0(jj),ml(jJ)/ml0(jj)
          enddo

         is=i
          call calhnew (ta,ml)
            call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
c              print*, 'repeat check', xkelvin*aw, awshell(I)*xkelvin,rh



           if (dabs( xkelvin*aw- rh).ge. 2D-3) then
           if (Idiagnose.ge.1) print*, xkelvin*aw, rh,i

c

           call aw_back
     & (ta,ML0,aw0,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           if (Idiagnose.ge.1) print*, 'before CALL_ML', aw0* xkelvin,rh


         call calhnew (ta,ml)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           
        if (Idiagnose.ge.1)     print*, 'after CAL_ML', aw* xkelvin,time


                    DO kk=1,np
                       ml(kk)=ml0(kk)
                       enddo

 
                       aw0=xkelvin*aw0
                       

           if (xkelvin*aw.lt.RH) then
           ff=1.0001
           else
              ff=1/1.001
           endif

            rha=rh/xkelvin
            DO Ii=1,90000
               DO kk=1,np
                  ml(kk)=ml(kk)/ff
                  enddo
               ml(1)=ml0(1)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

             if ( dabs(aw-rha).le..02D-3) goto 449
             if (aw0.ge.rh .and.aw.le.rha) goto 449
             if (aw0.le.rh .and.aw.ge.rha) goto 449
               enddo
 
 449           continue
        xn(I,1) = (xn(I,2)+xn(I,16)+xn(I,3)+xn(i,4)+xn(i,5) )*ML(1)/
     &(ML(2)+ml(16)+ML(3)+ml(4)+ml(5))

        awshell(I)=aw
c
        if (Idiagnose.ge.1)print*, 'finish again',aw*xkelvin,rh
C        print*, 'finish again',aw*xkelvin,rh
C        print*, 'finish again',aw*xkelvin,rh
c          
        if ( dabs(aw-rha).GE..01) THEN
           print*,'error: aw did not converge,stop'
           do J=1,NP
              WRITE(84,*) j, ML0(j)
           ENDDO
           STOP
        ENDIF
        
        
              endif


           
        endif
                endif

             enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      

            dtimesplitmax=time*.01 
       if (dtimesplitmax.le.1E-4) dtimesplitmax=1E-4
       if (dtimesplitmax.le.1D-3  .and. time.ge..05) dtimesplitmax=1E-3
ccccccccccccccccccccccc

       
         if (dtimesplitmax.ge. 1d0) dtimesplitmax=1d0



               if (dtimesplit.ge.dtimesplitmax) then
         if (imode_shell.eq.0) then
       call reset_shells(NS,x)
          endif
       endif
       
          xn22=0d0
               DO I=1,ns
                  xn22=xn22+xn(I,2)
                  if (x(I+1)-x(i).le. 0.6*dxmin)dtimesplitmax=-1
               enddo
c               print*,'dtimesplit',dtimesplit,dtimesplitmax
c            print*,'xn22, time' , time, xn22   
c               print*,'dtimesplit',dtimesplit,dtimesplitmax
               if (dtimesplit.ge.dtimesplitmax) then
       if (imode_shell.eq.1) then
c       call reset_shells(NS,x)


c           print*, 'imdoe_MA', imode_NH4NO3,imode_MA
         if (  imode_NH4NO3+imode_MA.ge.1 ) then

            dx = x(NS+1)-x(NS) 
            dx1 = x(NS+1)-x(NS) 
c            print*,'add b', dx,dxmin
            if (NS.ge.2) then
            dxx= x(NS)-x(NS-1) 
            if (dxx.ge.dx) dx=dxx
         endif
            if (NS.ge.3) then
            dxx= x(NS-1)-x(NS-2) 
            if (dxx.ge.dx) dx=dxx
         endif
         
            dd= dxmin*2
            
c     if (dd.le.10D-7) dd= 10D-7
c          print*,'ns00',NS
c            print*,'dddd',dx,dd
            rr=xn(1,2)/(xn(1,18)+ xn(1,2)+xn(1,29))
            rhh = .1 +0.4* rr
            
c            write(6,'(A,7E15.6)') 'rr',time, rr,rhh,rh
            
            if (RH.le. rhh  .and. dx1 .ge. dd)then
               if (idiagnose.ge.1)print*,'before addbin',time,ns
            call addbin(xn,x,NS,dd) ! split or merge bins
             if (idiagnose.ge.1) print*,'after addbin',ns
c          call split_shells_NH4NO3(x,ns) !NH4/NO3  gradient

          call split_shells_NH4NO3(x,ns) !NH4/NO3  gradient
          call  split_factor(x,NS)

         endif

          call merge_NH4NO3(x,NS)

c          print*,'ns0',NS
c     redistributue shells if uniform
c          print*,'nsA',NS
          if (rh.le.rhh .and.ilate.eq.0) ilate=1
          if (rh.ge.rhh+0.05.and.Ilate.eq.1)
     * call reset_shells_NH4NO3(NS,x)

          if (rh.ge.rhh+.1 .or.dl_ion.ge.1D-10)
     &          call reset_shells_NH4NO3(NS,x)

          
          if (idiagnose.ge.1) print*,'ilate', ilate, rh
          
       endif

      if (imode_NH4NO3+imode_MA.eq.0)   call merge(x,NS)
      if (imode_NH4NO3+imode_MA.eq.0)   call merge_redis(x,NS)
      
         call split_shells_aw(x,ns) !aw gradient

         if (imode_pH.le.1 .and. imode_MA.eq.0) call split_shells(x,ns) !pH gradient


      endif
 

      dtimesplit=0d0
            call cal_flux(time,x,flux,dtime,NS)

      endif
           
       

      
c     xloop make a output for 250000 timestep
      loopmax=1001
      if (idiagnose.ge.2) loopmax=1
      
        if (loop.ge.loopmax) then
c     diagnostic output finds the time limiting shell and species

           xmm=-1
           NN=Nmin
           if(NN.gt.NS) NN=NS
           if ( Jmin.le.np ) xmm= xn(nn,jmin)*55.1/xn(nn,1)
        if ( Jmin.ge.2013) xmm= xn(nn,jmin-2000)*55.1/xn(nn,1)
      if(jmin.eq.1001)
     &xmm= dlog(ml6shell(Nmin)*gamma2(Nmin,6))/dlog(.1d0)

      phh= dlog(ml6shell(Nn)*gamma2(Nn,6))/dlog(.1d0)
      IS=NS

      xcc=0
      DO j=1,np
         ml(J)=ml(1)*xn(NS,J)/xn(NS,1)
      enddo
      
       call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)

       DO j=1,np
c         ml(J)=ml(1)*xn(NS,J)/xn(NS,1)
         xcc=xcc+ izc(j)*ml(j)

      enddo
      xnp=0
      DO I=1,ns

         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         

      enddo
      if(idiagnose.ge.1) then
      write(6,'(4I7,10ES15.6)')nmin, jmin, isolidmin,loopkka,
     & time,sdtime/loop,xmm,phh, x(nmin+1)-x(nmin),xnp
         endif
         nmin1=nmin
         if (nmin1.ge.100) nmin1=nmin1-100
         write(26,'(4I7,10ES15.6)')nmin, jmin, isolidmin,ishelleq,
     &        time,sdtime/loop,xmm,phh, x(nmin1+1)-x(nmin1),
     & xn(Nmin1,6)/xn(nmin1,1)*ml(1), xn(Nmin1,18)/xn(nmin1,1)*ml(1)
         nn=NS-3
         if (NN.le.1) nn=1
        write(81,'(16E13.4)') time, (x(I+1)-x(I),i=NS,NN,-1),
     &  ((xn(I,12)+xn(I,29))/xn(I,1)*55.51,i=NS,NN,-1),       
     &        (xn(I,2)/xn(I,1)*55.51,i=NS,NN,-1)       ,rh,
     &        gamma2(ns,12)*gamma2(ns,18)

c        write(82,'(16E15.6)') time,   gamma2(ns,12),gamma2(ns,18)
        
     
      call flush(26)
      loopkka=0
      
      call flush(6)
      
      ss=0

          DO I=1,NS
c       phh = dlog(xn(I,6)*ml(1)/xn(I,1)*gamma2(I,6))/dlog(.1d0)
c      phshell(I)=phh


      DO J=1,NP
         ML(J)=ml(1)*xn(I,j)/xn(I,1)
      enddo
              IS=I
           call calHNew(Tdrop,mL)
          call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa1)
       phh = dlog(ml(6)*gamma2(I,6))/dlog(.1d0)
      phshellI=phh
      ss=ss+phh

c     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
 
      write(45, '(1I5,15E15.6)') I,time, awshell(I), phshelli
     & ,x(I+1)-x(I)
     &, xn(I,12)*ml(1)/xn(I,1)  
     &     , xn(I,2)*ml(1)/xn(I,1)    ,
     & (xn(I,29))*ml(1)/xn(I,1),dl_ions(I)
     & , xn(I,1), flux(I+1,1)-flux(I,1)
      
      if (awshell(I).ge.0d0) goto 223
      if (awshell(I).le.0d0) goto 223
      print*,'error: aw NaN stop'
      stop

 223   continue
      enddo
      call flush(45)
c      write(46,*) time, ss/NS
      call flush(46)
      sdtime=0
           loop=0
              endif


       sstime=sstime+dtime
       svap4=svap4+vap(4)*dtime
       svap18=svap18+ vap(18)*dtime
       sml12= sml12+xn(NS,12)/xn(ns,1)*mL(1)*dtime
       sml2= sml2+xn(NS,2)/xn(ns,1)*mL(1)*dtime
       sml18= sml18+xn(NS,18)/xn(ns,1)*mL(1)*dtime
       

      
      
      if (time.ge.outputtime(Nout)-1D-30) then
c         print*, 'dliq', dliq, xloop

c     in order to recalculate x for output
          call cal_flux(time,x,flux,dtime,NS)
         
         xloop=0
         
         if (xloop.ge.250000d0)then
            xloop=0
            DO k= Noutput,nout,-1
               outputtime(k+1)=outputtime(k)
            enddo
            
            Noutput=Noutput+1
            nout=nout+1
            print*, 'n',noutput, outputtime(nout)
            
         endif
       xhs=0
       
c     write titer concentration
       xn4=0d0
       xn18=0d0
       DO I=1,ns
          xn4=xn4+xn(I,4)
          xn18=xn18+xn(I,18)
       enddo
       DO J=1,np
       xn1(j)= 0
       DO I=1,ns
          xn1(j)=xn1(j)+xn(i,j)
       enddo
      enddo
      
      xnAN = (xn1(12)+xn1(18))+xn1(29)+xn1(30)+xn1(31)
      
      ratio = xnan/(xn1(2)+ xnan)
      
      DO J= 2,np-npsolid
         ml(j)= ml(1)* xn1(j)/xn1(1)
      enddo
      
       if (xn40.lt.0d0) then
          xn40=xn4
          xn180=xn18
       endif



      call aw_back(tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        

      xmv=xvion
              call cal_dl_an_suc(Tdrop,ML,aw,xmv,dli)           
      xmv=mv(29)
              call cal_dl_an_suc(Tdrop,ML,aw,xmv,dl29)           
              
      write(36,'(16E15.6)') time, aw, ratio, x(NS+1), dli,dl29,tdrop,rh
     &             , xvion
      
       
       
       if (time.ge.1D3) then
        write(3,'(2E15.6,I6,1E16.6,I5)') time,x(NS+1),ns,x(1),NP

          else
        write(3,'(1F17.7,1E15.6,I6,1E16.6,I5)') time,x(NS+1),ns,x(1),NP
        endif

        xis=0
        xns1=0
        DO kk=np-npsolid+1,np
           xns1=xns1+xn(1,kk)
           xis=xis+ mv(kk)*xn(1,kk)
           enddo
           xis=(xis/4/pi*3)**(1/3d0)




          write(6,'(I10,1E15.6,1E15.6,1F10.4,1E15.7)') 
     *         Nout,TIME,dtime,rh, awshell(NS)
          
          

          
             DO J=1,NP
             xns(J)=0
             enddo
             sna=0
             sco2=0
             sx=0
             scl=0
             sy=0
             xn11=0
             sflu=0
             ssars=0
             
             xcharge=0
             
             DO I=1,NS

                sna=sna+xn(I,16)
                scl=scl+xn(I,17)
                sco2=sco2+xn(I,5)
                sx=sx+xn(I,19)
                sy=sy+xn(I,20)


                xn11=xn11+xn(I,11)
                                sflu =sflu +xn(I,21)
                                ssars=ssars+xn(I,22)

                ML(1)=1000d0/18d0


             DO J=1,NP
                ML(J)=xn(I,j)* ML(1)/xn(I,1)
                xns(J)=xns(J)+ xn(I,J)
             enddo
c                      call ml2m(ML)



             DO kk=1,np
                      M(kk)=ML(kk)
                      enddo
                      is =I
       call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,POA,pl)
                      is =I
       call calHNew(Tdrop,mL)
      call aw_back(tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
       Smisch=ML(16)*ML(17)*gammaNa*gammaCl/Apnacl(tdrop)
       phh = dlog(ML(6)*gammah)/dlog(.1d0)
       
       DO J=6,NP
          xn(I,j) = ML(J)/ML(1)*xn(I,1)
            if(j.ge.8)    xcharge=xcharge+izc(j)*ml(j)
       enddo
       
       ppb=pacetic/press * 1D9
       vv=4*pi/3 * (x(I+1)**3-x(I)**3)


              flu=xn(i,21)/vv/ctiter0

              if (flu.le.1D-20) flu=1D-20
              sarsa=xn(i,22)/vv/ctiter0
              if (sarsa.le.1D-20) sarsa=1D-20
              
              


              DO KK=1,NP
                 M(kk)=ML(kk)
                 enddo
              if (m(18).le.1D-20) m(18)=1D-30


              
              xis=x(I)**3*4*pi/3


              Do kk=np-npsolid+1, np
              xis=xis + mv(kk)*xn(I,kk)
              enddo
             xis= (xis/4d0/pi*3)**(1/3d0)

!     radius fo shell I 

c     X(I+1) is the size with solid and liquid of shell I
c     xis is used for plot program
c     
c     output molalites,aw, radius and pH etc. in each shell 
c             phh=dlog(xhshell(i)* ml6shell(i))/dlog(.1d0)
                    DO KL=1,NP
                       M4(kl)=M(kl)
                       enddo
                       DO j=1,NP
c     idl has a problem to read tiny numbers
c                if (M(J).le.1D-30 .and.m(j).gt.0d0)  m4(J)=1D-30
             enddo
             DO k=1,npsolid
             s4(kk)=sshell(i,kk)
c             if (s4(kk).le.1D-30 .and. s4(kk).gt.0) s4(kk)=1D-30
             enddo

             if (I.eq.nS) then
       m4(12)= sml12/sstime
       m4( 2)= sml2/sstime
       m4(18)= sml18/sstime

                
             endif

             write(3,'(I6,501E14.5)') I,
     &        x(I+1),aw,(M4(j),J=1,20),flu,sarsa,(M4(j),J=23,NP),
     &   xhshell(i),pHH,xis,(sshell(I,kk),kk=1,npsolid), x(I+1)-xis
             
          




             
c     output the moles in each shells 
       if (time.ge.1D3) then
          write(23,'(I5,1E15.6,501E19.11)') i,time,
     &        (xn(i,j),J=1,np)
          write(21,'(I5,1E15.6,501E19.11)') i,awshell(i),
     &        dl_ions(i),time

       else
          write(23,'(I5,501E19.11)') i,time,
     &        (xn(i,j),J=1,np)
          write(21,'(I5,501E19.11)') i,awshell(i),
     &        dl_ions(i),time
          endif

          
c       write(43,'(I5,1E15.6,501E19.11)') i,time,gamma2(I,6)
c     &      ,gamma2(I,24), gamma2(I,25), gamma2(I,26)
       enddo
       
       


        if (imode_output.eq.4 .or. imode_output.eq.3) 
     & print*, 'phh  =' , phh

c     Titer convert to relative concentration
      sflu= sflu/stiter0
      if (sflu.le.1D-99) sflu=1D-99
      ssars= ssars/stiter0
      if (ssars.le.1D-99) ssars=1D-99
      
      xnna=0
      xnsolid=0
      xcharge=0
      xnoa=0
      xn22=0
      xnp=0
      xn55=0
      xn18=0
      xnca=0
      do I=1,NS
         NPL=np-npsolid
         xnna=xnna+ xn(I,16)+xn(I,npl+3)+xn(I,npl+4)*2
     & +xn(I,npl+6)*2 +xn(I,npl+9)+xn(I,npl+8)
         xn55=xn55+xn(I,5)
         xn22=xn22+xn(I,2)
         xn18=xn18+xn(I,18)
         xnoa=xnoa+ xn(I,29)+ xn(I,30)+ xn(I,31)
         xnp=xnp+ xn(I,24)+ xn(I,25)+ xn(I,26)+xn(I,npl+6)+xn(I,npl+8)
     &+xn(I,npl+11)         
         xnca=xnca+ xn(I,36)+ xn(I,npl+10)+xn(I,npl+11)         
      DO K=np-npsolid+1,np-npsolid+5
      xnOA= xnOA+xn(I,k)
      enddo
         
         
         DO K=np-npsolid+1,np
      xnsolid= xnsolid+xn(I,k)
      enddo
      
      enddo
c     add gas phase oxalic acid
      if (isOAclose.eq.1) then
        xngas= partvapoa*100/xn_aero /8.314/T 
        xnoa=xnoa+ xngas
      endif

      xn1=0d0
      DO I=1,NS
         do J=1,NP
            xn1(J)=xn1(j)+ xn(i,j)
         enddo
      enddo
      m(1)=1000/mm(1)
c     
      DO JJ=2,np
            m(JJ) = xn1(JJ)/xn1(1)*m(1)
         enddo

       call aw_back(tdrop,M,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
        call cal_MV(aw,tdrop,mv)
         
            rhoo= rho_mix(t,m,mm,mv)
            xl=0.59
       wAN =  m(12)*mm(12) +m(18)*mm(18)

       wt_suc = m(2)*mm(2)/(1000+m(2)*mm(2)+wAN)
       wt_AN =  wAN/(1000+m(2)*mm(2) +wAN)
             
        Ri=RI_ANSuc(wt_suc,wt_AN,xl,T,rhoo)
        ww=wAN/(wAN+ m(2)*mm(2) )
        alpha= alpha0+ (alpha1-alpha0)* ww
        
        write(34,'(7E15.6)') time,ri,rhoo,alpha
        
        
c        if (time.ge. timetr0(1)+ 100) stop
        
        
      if (time.ge. 1D3) then
      write(33,'(24E16.8,14F15.5)')
     &     time,sflu,ssars,xnsolid,tdrop,t,xnna,XNOA,xnp,xn18,xn22
         else
      write(33,'(1F17.7,113E16.8,14F15.5)')
     &     time,sflu,ssars,xnsolid,tdrop,t,xnna,XNOA,xnp,xn18,xn22
      endif
c      write(79,*) time, xn(1,24)+xn(1,25)+xn(1,26)+xn(1,37)
       IS=NS
       DO J=1,Np         
         ml(J)=xn(ns,j)/xn(ns,1)* ml(1)
        enddo
        call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,POA,pl)

        vap(4) =svap4/sstime
        vap(18) =svap18/sstime
        xno=0d0
        xnin=0d0
        xnnak=0d0
        DO I=1,NS
           xno=xno+xn(I,2)+xn(I,11)

           DO k=6,np-npsolid
              if (izc(k).ne.0d0) xnin=xnin+ xn(I,k)
              if (k.eq.16 .or. k.eq.19) xnnak=xnnak+ xn(I,k)*2
           enddo
           ratioo= xnin /(xnin+xno)
           ratioNa= xnnak /(xnnak+xno)
        enddo
c     output the vapor pressures
      if (time.ge. 1D3) then
         WRITE(14,'(15E15.6)') TIME, aw, vap(3)/press*1D9
     &           ,1D9*vap(4)/press,vap(18)*1D9/press,vap(17)*1D9/press,
     &         vap(5)/press*1D6
     &     ,tdrop,vap(29)/press*1D9,vap(33)/press*1D9,ratioo,rationa


c     &     ,1D9*pnh3/press,pHNO3*1D9/press,PHCL*1D9/press,PCO2/press*1D6
c     &     ,tdrop,pOA/press*1D9,pl/press*1D9
         else
         WRITE(14,'(1F17.7,14E15.6)') TIME, aw, vap(3)/press*1D9
     &           ,1D9*vap(4)/press,vap(18)*1D9/press,vap(17)*1D9/press,
     &         vap(5)/press*1D6
     &     ,tdrop,vap(29)/press*1D9,vap(33)/press*1D9,ratioo,rationa
      endif
      
c     mass solutes without crystal
       sstime=0d0
       svap4=0d0
       svap18=0d0
       sml12= 0d0
       sml2= 0d0
       sml18= 0d0
c
          
c     xms : mass of all solutes inclusive solid
c     xmw: mass of H2O
      xmw=0
      xmsl=0
      DO I=1,NS
         xmw=xmw + mm(1)* xn(I,1)
         xmsl=xmsl+xn(I,2)*mm(2)
         DO J=8,np-npsolid
         xmsl=xmsl+xn(I,j)*mm(j)
         enddo
         
      enddo
      xms=xmsl

      DO I=1,NS
         DO J=np-npsolid+1,np
         xms=xms+xn(I,j)*mm(j)
         enddo
      enddo
      
          rhoo=(xms+xmw)/4/pi* 3 / x(NS+1)**3

c     output the partial pressures
      if (time.ge. 1D3) then
          WRITE(15,'(15E15.6)') TIME, rh, partvap3/press*1D9
     &         ,1D9*partvap4/press,partvapHNO3*1D9/press,
     &         PartvapHcl*1D9/press,PartvapCO2/press*1D6, xms/(xms+xmw)
     & ,xmsl/(xmw+xmsL),awin, PartvapOA/press*1D9, Partvaplac/press*1D9

     
       else
          WRITE(15,'(1F17.7,12E15.6)') TIME, rh, partvap3/press*1D9
     &         ,1D9*partvap4/press,partvapHNO3*1D9/press,
     &         PartvapHcl*1D9/press,PartvapCO2/press*1D6, xms/(xms+xmw)
     & ,xmsl/(xmw+xmsL),awin, PartvapOA/press*1D9, Partvaplac/press*1D9


       endif
       

c     write projekted area and projected radius im cm2 and cm
       rr = x(1) + (x(NS+1)-x(1))/2d0
       write(18,'(3E15.6)') time,  rr*dsqrt(2d0),2*pi*rr**2 
       
             
           if (dliq.lt..1D-7) then

              print*, 'Only < 0.1 nm liquid coating left. '
              print*, 'It is then consider as fully solidfied, stop  '
              stop
           endif
      call flush(3)
      call flush(35)
      call flush(33)
      call flush(34)
      call flush(14)
      call flush(15)
      call flush(23)
      if (ssars.le.1D-99  .and. imode_output.ne.2 )  then
         print*,' finito, all the virus were inactivated stop'
         stop
      endif
      
                    Nout=nout+1

                    
          endif
c     end output
         if (idiagnose.ge.4) print*,'aaa7'
         if (idiagnose.ge.4) then
            DO J=1, npl 
               write(6,'(I5,100E15.6)') J,
     & (xn(I,j)/xn(I,1)*ml(1), i=1,ns)               
            Enddo

         endif
          

       DO KS= 1, NPsolid
           
           
          xnii=0d0
          DO I=1,ns
             xnii=xnii+ xn(i,np+ks-npsolid)
          enddo
          
          if (KS.eq.npsolid) then
c                      print*,'xnii', time, xnii
                      
                   endif

          if ( xnii.le.1D-60) then


             xapmax=0
             xapmin=1D10
             ns1=ns
             if (NS.ge.3) NS1=NS-1
             DO JJ = 1,NS1


                  if (xapmin.ge. sshell(jj,ks))xapmin=sshell(jj,ks)

                  if (xapmax.le. sshell(jj,ks)) then
                  imax=jj
                  xapmax=sshell(jj,ks)
                  xapmax1=sshell(jj,ks)
                  if (JJ.ge.2) xapmax1 =sshell(jj-1,ks) 

               endif
               
                enddo
                imax1=Imax-1
                if (imax1.lt.1) Imax1=1
                          DO jj=1,ns
                             xnsolidshell(jj)=0d0
                             DO J=np-npsolid+1,np
                             xnsolidshell(jj)=xnsolidshell(jj)+xn(jj,j)
                                enddo
                          enddo
                
                          IC=0
                          xsatic=0d0
                          NS1=NS
                          if (NS.ge.3) NS1=NS-1
                          DO jj=1,NS1
                           if (xnsolidshell(jj).gt.1D-30) then
                              
               if (xsatic.le. sshell(jj,ks)) xsatic=sshell(jj,ks)
                              IC=IC+1
                              ishell(IC)= I

c                              else
c                              xsatic= xapmax
                              endif
                         enddo

                         if (Ic.le.0d0)  xsatic= xapmax
                
                         isnuc=0
                         
                seff= AP_eff(KS)

c                if (xapmax.gt. 1. .and. xnsolidt.ge.1D-29 ) then
c     heterogeneous nucleation, half supersaturation
c                     seff =dsqrt(seff)
c                   endif



c     require shells with solid saturated
c     and the layer below is saturated
                   if ( xapmax .ge. seff  .and. xsatic.ge.1 .and.
     &               xapmax1.ge.1) isnuc=1



c     only for NaCl
                   if (inuc.lt.1) inuc=1
                   
                   timenuc_drop=timeeff(Inuc)


c                   if (ks.eq.9) print*,'inuc', inuc, timeeff(inuc)
c     NaCL and KCl nucleation at given times
                   
c             if ( (KS.eq. 9 .or. ks.eq.7).and. timenuc_drop.gt.0d0) then
             if (  timenuc_drop.gt.0d0) then
                      
                  if( timenuc_drop.ge.1 )   isnuc=0
               if (time.ge.timenuc_drop .and. timenuc_drop.gt.0d0) then
                  isnuc=0
                  print*, 'timenuc',  timenuc_drop, xapmax
                  
                  if ( xapmax.ge.1.01) then
                     isnuc=1
                      inuc=inuc+1
                      print*, 'isnuc', isnuc
                   endif
                   
c                  stop
                      
                   endif
                   endif

                   if (KS.eq.npsolid) then
c                      print*,'tt', time, inuc, xapmax
                   endif


                  
                   if (isnuc.eq.1) then
                    ff=1d0
                 if( timenuc_drop.ge.0 .and. KS.eq.9) ff=solidm(Inuc-1)
                   write(6,'(A,1I5,3E15.6)' ) 'solidm', KS, ff

                   
                      print*,'imax', imax, xapmax,xapmax1
                          xnca=0d0
                          DO II=1,NS
       xnca=xnca+ xn(iI,36)+ xn(Ii,npliq+10)
     &+xn(Ii,npliq+11)         
                       enddo
                       print*, 'befor P', xnca
                       
            xnsolido=0d0
            
            ML(1)= 1000/MM(1)


            DO Ii =1,NS
               DO K=2,NP
                  ml(k)=ml(1)* xn(ii,k)/xn(ii,1)
               enddo

       DO KK=1,NP
               M(KK)=ML(KK) 
       enddo
               IS=ii
               isaho = KS
               
               call  effl_ao(Tdrop,mL,msalt)               
               print*, 'msalt OA', ii,ks,msalt,sshell(ii,KS)
               

               
                  xnsolidJ = 0d0
              xjshell(ii)=0d0
               
               if (msalt.gt.0d0) then
              xnsolidJ = Msalt/ML(1)*xn(ii,1)*ff
              xnsolido=xnsolido+xnsolidj
              xjshell(ii)=xnsolidj

c     subtract component
              
              DO K=1, nsp(KS)
                 j= nsp_index(ks,k)
                 xnv=xnsp_nv(ks,k)
                 print*, 'j',j, xnv,xnsolidj
                 xn(ii,j) = xn(ii,j) - xnv*xnsolidj
              enddo
              
              
              DO kk=1,NP
              ML(kk) = ML(1) *xn(ii,kk)/xn(ii,1)
           enddo

              
              
              IS=Ii

              call calHNew(Tdrop,mL)
           call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           
           print*, ' saturation ',ii, ks, ss_solid(ks,tdrop,ml)
           print*, ' pHh solid ', dlog(ml(6)* gammah)/dlog(.1d0)
           
           
           endif
                    enddo
c                    stop
               print*,time, 'Solid efforesed !',KS,xnsolido

               
               if (xnsolido.le.0) goto 678
c

                   xnsolido= xnsolido
                   
                       

                       if(xnsolidt.gt.1D-30) then
                     
c     distribute xnsolido on the NaCl shells
                          DO I=1,ns
                             xnsolidshell(i)=0d0
                             DO J=np-npsolid+1,np
                             xnsolidshell(i)=xnsolidshell(i)+xn(I,j)
                                enddo
                          enddo
                          IC=0
                          
                        DO I=1,NS
                           if (xnsolidshell(i).gt.1D-30) then
                               IC=IC+1
                               ishell(IC)= I
                            endif
                            enddo
                           print*,'IC' ,IC
                           
                           ss=0

                           if (Ic.ge.1) then

                           Ishell(IC+1)=NS+1
                           
                           DO I=1,IC
                              I1= ishell(I)
                              if (I.eq.1) I1= 1
                              
                              I2=ishell(I+1)-1
                              ss1=0d0
                              DO J=I1,i2
                              ss1=ss1+xjshell(J)
                              enddo
                              print*, I,I1,i2, ss1,xn(I1,NP-npsolid+ks)
     & ,NP-npsolid+ks
                              xn(ishell(I),NP-npsolid+ks)=ss1
                           enddo
                           
c                              if (xnsolidshell(i).gt.1D-30)
c     & = xnsolidO/IC



                           goto 329
                    
                            endif

                         endif  !end xnsolid
                         

                       if(iscenter.eq.1) then
                          xn(1,np-npsolid+KS) = xnsolido
                          goto 329
                       endif

             if(iscenter.eq.0) then
c     distribute the solid into shells
c     with a spacing of xdissolid
                
                xx=0d0
                
c     the first shell
           JJ=1
                ishell(JJ)=1

c     out the solid only in shell with s > 1
                
                I=2

                NSS= x(NS+1)/xdissolid
                if (nss.le.1) NSs=1
                print*,'nss', nss,xnsolido

                if (Nss.eq.1) then
                  
                   jj=Ns
                   if (jj.ge.3) jj=ns-1
                    
                   xn(jj,np-npsolid+ks) = xnsolido
                   print*,'jj', jj,ns
                   goto 329
                endif

                
           x11= x(1)+ 0.05d0*(x(2)-x(1))

c           x11= (x(1)**3+ xnsolido/Nss*mv(NP-1)/4/pi*3)**(1/3d0)

                ss=x11**2
 309            continue
                    dx =  (x(I+1)-x11)

                write(6,'(I5,5E15.6)')
     & I ,dx*1D4


c         if (dx .ge. xdissolid .and. sOAshell(I).ge. 0.8 ) then
         if (dx .ge. xdissolid ) then

                    JJ=JJ+1
                    print*,' jj ', jj,I
                  

                       x11=x(I) +.05d0*(x(I+1)-x(I))


                    ishell(JJ)=I
                   ss=ss+ x11**2
                                     endif
                   I=I+1
                   if (I.LE.NS) goto 309

                   
                   
             print*,'solid and total bins ', jj,NS
             if (Ishell(jj).eq.ns.and.jj.gt.1.and.ishell(jj-1).lt.NS-1)
     &            ishell(jj)= nS-1
             
             IC=JJ
                           Ishell(IC+1)=NS+1
                           ss=0
                           DO I=1,IC
                              I1= ishell(I)
                              I2=ishell(I+1)-1
                              ss1=0d0
                              DO J=I1,i2
                              ss1=ss1+xjshell(J)
                              enddo
                              print*, I,I1,i2, ss1
                              xn(I1,NP-npsolid+ks)=ss1
                              ss=ss+1
                           enddo
                           
                           print*, 'xnsolid ', xnsolido, ss
c                snn=0
c             DO I=1, JJ
c                II= ishell(I)
c                   x11=x(Ii) +.05d0*(x(Ii+1)-x(iI))
c                xn(II,NP-npsolid+ks) =xnsolido/ss*x11**2
c                snn=snn+1/ss*x11**2
c                enddo

c     distributes the molecules weighted with the area
c                print*, ' normalize factoer =1 ', snn
                         endif !iscneter=0
                       
                       
 329                   continue


        DO I=1,NS
           DO J=2,Np
              ml(J)=ml(1)*xn(I,j)/xn(I,1)
           enddo
           call calHNew(Tdrop,mL)
           call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           
c           print*, ' pHh AHO ', dlog(ml(6)* gammah)/dlog(.1d0)
c           ss= gammaH2PO4 * gammana * ml(16)*ml(25)/apNaH2po4(T)
           ss= ss_solid(KS,tdrop,ml)
c           print*, 'ss', I, ss
           if (time.ge.timeeq)then

             xkelvin = dexp( 2* sigma * MV(1) /(8.314E7*Tdrop*x(NS+1)) )
          rhh= rh/xkelvin
         dd=.3D-3
         if (x(NS+1).le.0.1) dd=3D-3
         

        if (dabs(Aw-rhh).ge.dd) then

           IS=I
           aw1=aw*xkelvin
           call cal_ml(Ta,rhh,ML)
           

           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           loop=loop+1
           if (loop.ge.500) then

           loop=1
           endif
         DO kk=1,np
           Ml0(kk)=ML(kk)
           enddo

        xn(I,1) = (xn(I,2)+xn(I,16)+xn(I,3)+xn(i,4)+xn(i,5) )*ML(1)/
     &(ML(2)+ml(16)+ML(3)+ml(4)+ml(5))


c      repeat check
          ff= 1000d0 /(MM(1)*xn(I,1))

          DO J=2,NP
             ML(J)= ff* xn(I,j)
          enddo
         is=i
          DO kk=1,10
          call calhnew (ta,ml)
            call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
c c                      print*, 'repeat check', xkelvin*aw, rh
            enddo

           if (dabs( xkelvin*aw- rh).ge. .3D-3) then
              print*, xkelvin*aw, rh, xkelvin



           call aw_back
     & (ta,ML0,aw0,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           print*, 'before CALL_ML', aw* xkelvin
           print*, ml(29)+ml(30)+ml(31)
           print*, ml0(29)+ml0(30)+ml0(31)

         call calhnew (ta,ml0)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           print*, 'after CAL_ML', aw* xkelvin,time


                    DO kk=1,np
                       ml(kk)=ml0(kk)
                       enddo

 
                       aw0=xkelvin*aw0
                       

           if (xkelvin*aw.lt.RH) then
           ff=1.0001
           else
              ff=1/1.001
           endif

            rha=rh/xkelvin
            DO Ii=1,50000
               DO kk=1,npl
                  ml(kk)=ml(kk)/ff
                  enddo
               ml(1)=ml0(1)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
c             print*, ii, aw
             if ( dabs(aw-rha).le..02D-3) goto 484
             if (aw0.ge.rh .and.aw.le.rha) goto 484
             if (aw0.le.rh .and.aw.ge.rha) goto 484
               enddo
 
484           continue
        xn(I,1) = (xn(I,2)+xn(I,16)+xn(I,3)+xn(i,4)+xn(i,5) )*ML(1)/
     &(ML(2)+ml(16)+ML(3)+ml(4)+ml(5))
              
c
        print*, 'aw0,rh', aw0,rh
        print*, 'finish again',aw*xkelvin,rh,ii,ff
C        print*, 'finish again',aw*xkelvin,rh
C        print*, 'finish again',aw*xkelvin,rh

        if ( dabs(aw-rha).GE..01) THEN
           print*,'error: aw did not converge,stop'
           do J=1,NP
              WRITE(84,*) j, ML0(j)
           ENDDO
           STOP
        ENDIF
        
        
              endif


           
        endif
          
              
           endif
           
           
        enddo
        
        
                          xnca=0d0
                          DO II=1,NS
                          xnca=xnca+ xn(iI,36)+ xn(Ii,npliq+10)
     &+xn(Ii,npliq+11)         
                       enddo
                       print*,'after ca', xnca,ks

        
        endif !nuc
        

        
      endif !xii

      enddo

 678  continue
      
ccccccccccccccccccccccccccc
c     End efforencens oxalate
cccccccccccccccccccccccccc
               if (idiagnose.ge.4) print*,'aaa8'
         if (idiagnose.ge.4) then
            DO J=1, npl 
               write(6,'(I5,100E15.6)') J,
     & (xn(I,j)/xn(I,1)*ml(1), i=1,ns)               
            Enddo

         endif

      
c     goto time loop 11
          
          IF (TIME .LE. timetr0(ntr0) ) GOTO 11

          

 200    continue

        open(2,file='output_ML.dat')
c     output molalities of outermost shell
c     it serves for the input ML.dat

        print*,'finito'

c         xn1=0d0
           DO J=1, NP
           xn1(j)=0d0
        DO I=1,NS
           xn1(J)=xn1(J)+ xn(I,j)
        enddo
      enddo

      DO J=1,NP
      write(2,*) J,ml(1)*xn1(J)/xn1(1)
      enddo
      
            write(2,'(A,2E15.6)') 'aw = ', aw
            write(2,'(A,2E15.6)') 'pH = ', phh
            write(2,'(A,2E15.6)') 'pCO2 = ', pco2*1E3 !ppm
            write(2,'(A,2E15.6)') 'pAA = ', pacetic*1E6 !ppb
            write(2,'(A,2E15.6)') 'pNH3 = ', pNH3*1E6 !ppb
            write(2,'(A,2E15.6)') 'pHNO3 = ', pHNO3*1E6 !ppb
            write(2,'(A,2E15.6)') 'pHCl = ', pHCl*1E6 !ppb
            write(2,'(A,2E15.6)') 'pOA = ', poA*1E6 !ppb
            write(6,'(A,2E15.6)') 'RH = ', aw

             print*, 'endtime reached'



        close(2)

        

        
       stop
       end

      
c     ---------------------------------------------------
c     Functions and subroutines
c     ---------------------------------------------------


c     -------------------------------------------------
c
c     LOOKUP TABLE for water activity and activity coefffients

c      
c     -------------------------------------------------

           

c     -----------------------------------------------------------------------------------------------------
c     This subroutine calculate the flux f of all species and timestep dtime including solid NaCl in all shells.
c
c     x
c     f(i,np) is the flux of the inner liquid shell to the NaCl crystal
c     f(i,np+1) the flux from the liquid to NaCl solid in the shell i
c     f(i,j): is the flux of species j from shell i to i-1.
c     f(NS+1,NP) is the flux from gas phase to the particle
     
c     output x: NS+1 array, x(i) and x(i+1) the lower and higher boundary (radius) of first shell
c
c     -----------------------------------------------------------------------------------------------------
c     the time step is dynamically calculated to ensure stable compostion and pH
c      



c     ------------------------------------------------------------------------------------
c
c     calculate the initial number of moles of each species in all shells xn at given molalities and size r0
c     each shell has the same volume 
c
c     ------------------------------------------------------------------------------------

      subroutine set_shells_vol(t,RH,r0,ML,x,NS,NSmax)

       implicit real*8 (a-h,k,m, o-z)
       parameter (np=50,NSMM=100)
       real*8 MM(NP),MV(NP),dl_factor2(2,NP),ntiter(NP)
c       common /DL/ DL_factor2 ,deltaxgas
        common /DL/ DL_factor2,deltaxgas,jmin,nmin

       common /M/ MM,mv,izc(np)          
       real*8 x(*),xn(NSMM,np),ml(*)
       common /xn/xn
       common /rcore/ rcore

         xnn=(10D-7)**3/(0.2D-7)**3

         XNTITER0=xnn/6.023D23


 
         
       pi=dacos(-1d0)

        xv=ML(1)*MV(1)+ML(2)*Mv(2)
       DO I=6, np
          xv=xv+ ml(I)*mv(I)
       enddo
          xv=xv+ ml(NP)*mv(NP)
       xvmol =XV
       print*, xvmol

       

       
c     divide into shell with distances of about  0.3 - 2  nm.
c
c     calculates the diffusion length
         
       

c       dx = r0/Ns
       call cal_dlaw(T,rh,dl)
       dmin= dsqrt(dl*2)



        ns=nsmax


        if (NS.gt.NSmax) NS=NSmax
       if (NS.le.1) NS=1


       print*, NS

       
       v1 = 4d0*pi/3d0*r0**3/NS
       xntiter=xntiter0/ns
       
       if (ns.gt.NSMM) then
          print*, 'error: increase the dimension for the size bins '
          print*, ' or increase the minimum thickness!  '
          stop
          endif
        
c     Define the radius of the inner shells
c       DO I=2,ns+1
c           x(i)=dx*(I-1)
c           print*, I,x(I)
c        enddo

c     calculates the moles of H2O and organics in each shell
          vcore = 4*pi/3 * rcore**3
          x(1)= rcore
          DO I=1,NS
             ntiter(I)=xntiter
             
             v= v1*I
          x(I+1)= (vcore+v/4d0/pi*3d0)**(1d0/3d0)

          DO J=1,nP
          xn(I,j)= ML(J)*v1/xvmol
c          print*,j,mm(J),mv(J)
       enddo
          v2=xn(I,1)*mv(1)+xn(I,2)*mv(2)
          DO J=6,NP
             v2=v2+ MV(J)*XN(I,j)
          enddo
          print*,I, v1,v2
       enddo
       print*,  ntiter(NS)
       DO I=1,NS+1
          print*,i , x(I)
          
       enddo
       return
      end



c     ----------------------------------------------------------------------------------------
c
c     calculate the diffusion coefficient of ions as a function of aw
c
c     ----------------------------------------------------------------------------------------

      

c     aw: according to He (is not used)
c     dl: diffusion coefficient in cm2/s


      subroutine cal_dlaw(T,aw,dl)
      implicit real*8 (a-h,o-z)
      T0=293.15
      aw0=aw
c     take the t-Dependence of citric acid
            if (aw.ge.1d0) aw0=1d0
      call cal_ions(T0,aw0,dl)

c     takes the T dependence of citric acid
      call  cal_dlaw_citric(t,aw0,dlt)
      call  cal_dlaw_citric(t0,aw0,dlt0)
c      dlt0 = Dh2o_func(T0f
c      dlt  = Dh2o_func(T)
      dl= dl* dlt/dlt0

      return
      end


c     ----------------------------------------------------------------------------------------
c
c     calculate the diffusion coefficeint of H2O as a function of aw obtained from Figure 5 of Walker et al paper
c
c     ----------------------------------------------------------------------------------------

      subroutine cal_dlaw_walker(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)
      common /enhance/radius,r1,enh_factor
      
      common /Ienhance/Ienh,iscenter,isahoo
      

      tt=298.15
      dl0=dh2o_func(TT)

      enhance = 1


      ff= enhance
            fac(1)=7E-9*ff
            fac(2)=8E-9*ff   !10E-9
           fac(3)= 8D-9*ff !11D-9
           fac(4)= 9E-9*ff

       fac(5)= 1.8E-8
       fac(6)= 6E-8
       fac(7)= dl0/10
       fac(8)= dl0/2
       fac(9)= dl0

      rh(1)= 0
      rh(2)=.5
       rh(3)=.65
       rh(4)=.75

       rh(5)=0.85
       rh(6)=0.90
       rh(7)=0.96
       rh(8)=0.99
       
       rh(9)=1
       DO KK=1,9
      fac(kk)=dlog(fac(kk))
      enddo
      
      N1=1
       N5=9
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))
       t0=298
       call  cal_dlaw_citric(t,aw0,dlt)
       call  cal_dlaw_citric(t0,aw0,dlt0)
       dl= dl* dlt/dlt0
       
      

       return
      end



c     ----------------------------------------------------------------------------------------
c
c     calculate the diffusion coefficient of H2O as a function of aw obtained from from present EDB study used in present paper
c
c     ----------------------------------------------------------------------------------------
   
      subroutine cal_dlaw_walker_h2o(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)
      data fac1 /1d0/
      data fac2 /1d0/
      logical ex
      save fac1,fac2 , key,ex,fac3,fac4
      xv=18d0
       call  cal_dlaw_mod(T,aw0,dl,xv)

      return

      
      if (key.eq.0) then
      INQUIRE (FILE='DL_SLF.dat', exist=ex)
      fac1=1d0
      fac2=1d0
      fac3=1d0
      fac4=1d0
      key=1
      
      if ( ex ) then
         open(99, file='DL_SLF.dat')

         read(99,*,end=23) fac11
         read(99,*,end=23) fac11
         read(99,*,end=23) fac11
         read(99,*,end=23) fac11

         read(99,*) fac1
         read(99,*) fac2
         read(99,*,end=23) fac3
         read(99,*,end=23) fac4
         

 23      continue
         close(99)
      endif
      endif
      

      
      aw1=1

      t0=298.15
      Dl0=  2.44D-5  ! scale d0 to 2.44E-5 cm2/s at 298d0 and aw=1
      dl0= dh2o_func(T0)
      
            fac(1)=7E-9*5*fac1
            fac(2)=8E-9*5*fac2   !10E-9
           fac(3)= 8D-9*5*fac3 !11D-9

           fac(4)= 9E-9*5*fac4

       fac(5)= 1.8E-8*5
       fac(6)= 6E-8*5
       fac(7)= dl0/4
       fac(8)= dl0*.8
       fac(9)= dl0

      rh(1)= 0
      rh(2)=.5
       rh(3)=.65
       rh(4)=.75

       rh(5)=0.85
       rh(6)=0.90
       rh(7)=0.96
       rh(8)=0.99
       
       rh(9)=1
       DO kk=1,9
      fac(kk)=dlog(fac(kk))
      enddo
      N1=1
       N5=9
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))

c     temprature dependence
       
       t0=298.15
       call  cal_dlaw_citric(t,aw0,dlt)
       call  cal_dlaw_citric(t0,aw0,dlt0)
       dl= dl* dlt/dlt0

       return
      end

c     ----------------------------------------------------------------------------------------
c
c     calculate the diffusion coefficient of H2O of citric acid Liendhar et al
c
c     ----------------------------------------------------------------------------------------
      
      subroutine cal_dlaw_citric(t,aw,dl)
      implicit real*8 (a-h,o-z)

       a1=0.61477E+00
       a2=0.42200E+00
       a3=0.90000E-02
       b1=0.13800E+02
       b2=0.15693E+00
       b3= -0.90000E-02

c    7    -0.78000E+01
c    8   -0.20000E+02


         DH2O_1=10**(-6.514-387.4/(T-118.0))
         DH2O_0=10**(-15.0-175.0/(T-208.0))
         tc=t-273.15
         if (tc.ge.-7.8) tc=-7.8
         A=a1+a2*Tc+a3*Tc**2
         if (tc.ge.-20) tc=-20
         b=b1+b2*Tc+b3*Tc**2
         alpha=dexp((1-aw)**2*(A+aw*B))
!         print*, a,b
         dl=DH2O_1**(alpha*aw)*(DH2O_0**(1-aw*alpha)) ! m2/s
         dl=dl*1D4 ! cm2/s
         return
         end

      
c     ----------------------------------------------------------------------------------------
c
c     calculate the diffusion coefficient of H2O of sucrose Zobrist et al
c
c     ----------------------------------------------------------------------------------------

      subroutine cal_dlaw_suc(T,aw,dl)
         implicit real*8 (a-h,o-z)

       real*8 x(9)


      data x / 0.175,  -46.46,1.7,262.867, 10.53,-0.3,127.9,
     & 0.4514,-0.5/

c       omega1= 1-w

        

c      DO Aw= 0,1.01,.2
      x2= 1-aw
      if (x2.lt.0d0) x2=0d0


       if (x2.lt.0) x2=0

        a =  x(1)*(1+x(2)*x2)
        b =  x(4)*(1+x(5)*x2+x(6)*x2**2)
         t0 = x(7)*(1+x(8)*x2+x(9)*x2**x(3))
c         print*, aw,x2, t0
c         enddo





        dmin=1d-30

        if (t.le. t0+0.1d0) then
           dl=dmin
           else

       xx= a+ b/(T-t0)

       fcal = 10d0**(-xx)

       fcal=fcal*1d-7 ! m2/s

       dl =fcal
       endif
       if (dl.lt.dmin) dl=dmin
c      print*, aw, t0, dl
       dl=dl*1D4 !cm2/s
c      enddo

c
       return
       end




c     ----------------------------------------------------------------------------------------
c
c     calculate aw of sucrose: Zobrist et al.
c     
c     ----------------------------------------------------------------------------------------

       subroutine calaw_beni(T,omega1,aw)
       implicit real*8 (a-h,o-z)
       data a /-1/
       data b /-0.99721/
       data c /0.13599/
       data d /0.001688/
       data e /-0.005151/
       data f /0.009607/
       data g /-0.006142/
 

       w2=1-omega1
       T0 =298.15d0
       aw = (1+a*w2)/(1+b*w2+c*w2**2)
       aw =aw+ (T-t0)*(d*w2+ e*w2**2+f*w2**3+g*w2**4)

      
       return
       end







c     ------------------------------------------------------------------------
c   henry's law coefficient of 
c      H+(aq) + NH3(gas)    --> NH4+(aq)  
c      k =( a_H+ * pNH3 ) / a_NH4+
c     in bar-1 
c     Rennard 2004
      
            function xknh3(T)
      implicit real*8 (a-z)

      xkw =  -0.61205E+01+ 0.44820E+04/t+0.17055E-01*t
      xkw=10d0**(-xkw)
      xkb = dexp( 16.9732 - 4411.025/T -0.044*t)
      xH =  dexp( -8.09694 + 3917.507/T -0.00314*t)
c      print*, xkw, xkb, xh
c      Renard 2004


      xknh3 = Xh*xkb/xkw
      return
      end

c     ------------------------------------------------------------------------
c   dissociation  coefficient of 
c      H+(aq) + NH3(aq)    --> NH4+(aq)  
c      k =( a_H+ * aNH3 ) / a_NH4+
c     in M
c     Rennard 2004

      
            function xknh4b(T)
      implicit real*8 (a-z)

      xkw =  -0.61205E+01+ 0.44820E+04/t+0.17055E-01*t
      xkw=10d0**(-xkw)
      xkb = dexp( 16.9732 - 4411.025/T -0.044*t)
c      xH =  dexp( -8.09694 + 3917.507/T -0.00314*t)
c      print*, xkw, xkb, xh
c      Renard 2004


      xknh4b = xkb/xkw
      
      return
      end


c     --------------------------------------------
c
c     required by calHnew , calculated the H+ and OH- concentration to ensure to charge balance
c     
c     ------------------------------------------

      

      function fcnnew(x)

        IMPLICIT REAL*8 (A-H,O-Z)
       integer NP

       parameter (np=50,NSMM=100)
       real*8 ML(NP), NL(NP)
       real*8 mm(NP),mv(NP)
        integer izc(NP)
       common /M/ mm,mv,izc
       common/suls/ T, ML, ppartco2
c        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
c     &      ,gammahpo4 ,gammaHOA,gammaOA,gammah3po4
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak

        
        common /acid/ xkacid1, xkacid2,xhacid,xkacid3
        


       
       
       common /iskinetic/iskinetic

       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq
       
         real*8 gamma2(NSMM*2,NP)
        common /gamma2/gamma2
        data key /0/
        save key

        if  (key.eq.0) then
           key=1
           DO I=1,nsmm
              DO J=1,np
                 gamma2(i,j)=1d0
              enddo
                 enddo

                 
              endif


        
       
c       common /restart/ 

       HCO2 = 0.034 * dexp(2300 *(1/T-1/298.15))

      
c     new ML(23): molecular nh3: RG bates Pinching 1949
c     KNH4 (NH4+ + H2O - NH3 + H3O+
c     

c     NH3 + H+ --> NH4

c      xKNH4bb = 10d0**(9.4d0)
c       gammah
c       if (is.gt.NS .and. ns.ge.1) IS=NS

c     take index 1 

       if (is.lt.1 ) IS=1
       if (is.gt.2*NS) IS=1
       if (NS.le.0) IS=1
       

       
       I=IS
       gammaco2=gamma2(I,13)
       gammah=gamma2(I,6)
       gammaNH4=gamma2(I,12)
       gammaCl=gamma2(I,17)
       gammaS1=gamma2(I,27)
       gammaS2=gamma2(I,28)
        gammaoh= gamma2(I,7)

        gammaHCO3=gamma2(I,15)
       gammaCO3= gamma2(I,14)
       gammaH2PO4=gamma2(I,25)
       gammaH3PO4=gamma2(I,24)
       gammaHPO4=gamma2(I,26)
       gammaHOA=gamma2(I,30)
       gammaOA=gamma2(I,31)
        gammaAA = gamma2(I,18)! gammaNO3
        gammaLA = gamma2(I,18)! gammaNO3
        gammaPO4 = gamma2(I,37)! gammaNO3
        gammaA3 = gamma2(I,38)! gammaNO3


       
       aw = gamma2(I,1)

          xKNH4bb=xknh4b(T)



        ml(6)=x
        xsur=ml(2)
        XNH4= ml(4) !total NH4+ + NH4CH3COO
        xace= ML(3) !total acetic acid
        xco2=ML(5)
        xh= x
        xsul=0d0 ! no sulfat
        xno3= ML(18)

        
c     takes the activity coefficients of H+, NH4+ into account
c     

c

        fcn0=1D10
        itmax=50
        DO II=1,2

         xkw=  dexp( -0.92644d1-0.68727E+04/t) !dissociation of H2O
c         xkw=xkw/gammaH



c         xka=1.74D-5            ! dissociation constant of acetic acid !
c     Harned and Ehlers 1933

       xx= -1500.65/T-6.50923 * dlog(T)/dlog(10d0)-0.0076792* T
     & +18.67257d0
         xka= 10.d0**(xx)

        xka=xka/gammaH/gammaAA  ! takes the activity coefficient of H+

c     xksalt=10d0**(-9.53d0)    !  k for NH4CH3COO

        xksalt=4.E-5    !   Jaffe 1991 H =110 KJ/mol  and boiling point 117.1
c total =1.327E-2 hP
        xksalt=3.44E-3
        xksalt=0d0
        
        

        xksalt=xksalt/gammaNH4/gammaAA

        xkNH4bb=xknh4b(T)*gammaH/gammaNH4

c         av = xksalt *( 1 +1/(x*xkNH4bb))

c        bb=-xace+XNH4 + av*(1+x/xka) 
 
c        aminus = 1d0/2d0/(1d0+x/xka) *(-bb+
c     & dsqrt( 4d0 * (1d0+x/xka)*xace*av + bb**2d0))
        aminus= ML(3)/(1+ x/xka)
        HA= aminus*x /xka

        xsalz =0d0

c     dissociation NH3
        
c        XNH4plus = XNH4 - xsalz
        XNH4A = XNH4 - xsalz

        XNH4plus = XNH4A /(1 +  1/(x *xkNH4bb)) 

c     set to ML
        ML(10)=xsalz
        ML(9)= HA              !acetic acid
        ML(8) = aminus         !CH3COO-
        ML(12)= xNH4plus
        ML(23)= xnh4A-xNH4plus

        ml(6) =x
         
c       call aw_back(t,ML,aw,gammaHa,gammaNO3,gammaNH4,gammaCl,gammaNa)        
c       gammah= dsqrt(gammaH* gammaHA)

        aOH = aw *xkw/(x*gammaH)

c     takes the water activity of NaOH for OH-
        
        ML(7)= aOH/gammaOH
        
        
c        print*,'H+= ', x
c        print*,'HA= ', HA ! HA, A-
c        print*, 'A-  OH- ', aminus, xoh
c        print*, 'Asalt = ', xsalz, xnh4plus*aminus/xksalt
c        print*, 'NH4+= ',  xnh4plus
c        print*, 'xace, xnh4= ', xace, xnh4
c        Xk3 = 1.7D-3            ! CO2 + H2O --> H2CO3 at 298.15K
c        Xk4 = 2.5D-5            ! H2CO3 -> H+ HCO3- 

c     disscoci ation CO2(aq) + H2O --> H+ HCO3-

c        xk34 = 4.448E-7*dexp(-2133*(1/T-1/298.15)) ! 
      xk34= 290.9097 -14554.21/T-45.0575*log(T)
      xk34=dexp(xk34)

        xk34= 290.9097 -14554.21/T-45.0575*log(T)
      xk34=dexp(xk34)

      xk34=xk34/gammaH/gammaHCO3


        
          xm15 =xk34*aw*gammaco2/(x)

c      HCO3- --> H+ + CO3-2
c     reaction CO2(aq) + H2O --> H+ HCO3-  (1)
c     reaction CO2(aq) + OH- -->   HCO3-   (2)



c     disscociation HCO3- --> H+ + CO3-2

         xk2 = 10D0**(-10.33)*dexp(-3347.3*(1/t-1/298.15d0))
         xk2 = XK2*gammaHCO3/gammaH/gammaCO3
 
         xm14= xm15* xk2/x
        
        if (iskinetic.eq.0) then

        if (iseqco2.eq.0) then
           
        ML(13) = ML(5)/(1+xm15+xm14)   !CO2
        if (ml(13).le.0)ml(13)=0d0
        ML(15)= ML(13)*xm15! HCO3-
        ML(14) = ML(13)*xm14    !CO3-2
c           xss =ml(15)+ml(14)
c           ml(15) = xss/(1+ xk2/x)
c           ml(14)=xss- ml(15)

      else
c           print*,' ppartco2 ', ppartco2,hco2
           ML(13)=HCO2*ppartco2/1013.5d0

           ML(15) =xm15*ML(13)
           ML(14) =xm14*ML(13)
           if (ML(14)+ml(15).ge.30) then
              xmmm= ML(14)+ML(15)
              ML(14) = 30/xmmm * ml(14)
              ML(15) = 30/xmmm * ml(15)
                        endif


           ML(5)=ML(13)+ml(14)+ml(15) 
           endif

        endif
        
c           print*,'iseqco2', iseqco2
        if (iskinetic.eq.1) then
        if (iseqco2.eq.0) then

c        ML(13) = ML(5)/(1+xm15+xm14)   !CO2
c        if (ml(13).le.0)ml(13)=0d0
c        ML(15)= ML(13)*xm15! HCO3-
c        ML(14) = ML(13)*xm14    !CO3-2
           xss =ml(15)+ml(14)
           ml(15) = xss/(1+ xk2/x)
           ml(14)=xss- ml(15)
c           print*, ml(14),ml(15),ml(13)
        else
c           print*,' ppartco2 ', ppartco2,hco2
           ML(13)=HCO2*ppartco2/1013.5d0

           ML(15) =xm15*ML(13)
           ML(14) =xm14*ML(13)
           if (ML(14)+ml(15).ge.39) then
              xmmm= ML(14)+ML(15)
              ML(14) = 39/xmmm * ml(14)
              ML(15) = 39/xmmm * ml(15)

              endif


           ML(5)=ML(13)+ml(14)+ml(15) 
           endif
           endif

           


c           xkP1= 6.9E-3*gammah3po4/gammaH/gammaH2PO4

           xkp1 = 10d0**(799.31/T-4.5535+.013486*T) !bates 1951
           xkp1= 2.85829E4/T-501.5113 -0.094568*t+78.9595*dlog(t) !wright 1979
     &          - 18.654E5/T/T
c           print*,'k1', xkp1
           
           xkP1= dexp(xkp1)*gammah3po4/gammaH/gammaH2PO4


c           xkp2= 6.2E-8/gammaH/gammaHPO4*gammaH2PO4

           xkp2= 10**(-7.2)*( dexp(7100/8.314*(1/t-1/298.15))) !Vega Romero 1994
           xkp2=dexp(2.7434E4/T-485.0409-0.0938*T+75.0608*dlog(T)
     & - 20.6675E5/T/T)                
c           print*,'k2', xkp2
           xkp2= xkp2/gammaH/gammaHPO4*gammaH2PO4
c     
           
           xkp3= 10.d0**(-12.32d0)/gammaH/gammaPO4*gammaHPO4
           xkp3 = dexp(16482.1/T-858.75 -0.2883*t +151.09*dlog(T))

           xkp3= xkp3/gammaH/gammaPO4*gammaHPO4

           xm25 = xkp1/x
           xm26 = xkp2/x * xm25
           xm37 = xkp3/x * xm26

           xptot= ml(24)+ml(25)+ml(26)+ml(37)
           ML(24) = xptot/(1+xm25+xm26+xm37)
           mL(25) = ml(24)* xm25
           mL(26) = ml(24)* xm26
           mL(37) = ml(24)* xm37

           dh0=7036.629444310801
           xms = ml(27) + ml(28)
	xlk=-4.556380021818660+dh0*(1/298.15-1/T)-275./8.314*
     &  dlog(t/298.15d0)
	xks=exp(xlk)
        xks = xks /gammah/gammas2*gammas1
        

        ML(27) = xms /(1 + xks/ML(6))
        ML(28) = ML(27) *xks/ML(6)
c        print*,ml(27),ml(29),xms
        
        
c     Oxalix acid
        xmOA = ml(29)+ml(30)+ml(31)+ml(38)
        

c     Kettler 1998 table 1 Kurz + Farrar

         xka = 5.1E-2 
         xka2 = 5.E-5 
         fcal =  -0.13216E+01 -0.12564E+03*(1/t-1/298.15d0)

         xka= 10d0**fcal

c     kettler 1998 table 13
	fcal=-0.42640E+01 -0.31566E+04 *(1/t-1/298.15d0)
     & -0.11807E+02* dlog(t/298.15d0)
        xka2= 10d0**fcal
c        common /acid/ xkacid1, xkacid2,xhacid
        if (xkacid1.gt.0d0) then
        xka=xkacid1
        xka2=xkacid2
        xka3=xkacid3
        if(xka2.le.1D-50) xka2=0d0
        if(xka.le.1D-50)   xka=0d0
        if(xka3.le.1D-50)   xka3=0d0
        
      endif


        xka=xka/gammaH/gammaHOA  ! take the activity coefficient of H+
        xka2=xka2/gammaH/gammaOA*gammaHOA  ! take the activity coefficient of OA-1 OA-2
c        write(6,'(A,16E15.6)') 'AAAA', x, ml(7),aw,aoh,gammaoh,t,xkw
        
        xka3=xka3/gammaH/gammaA3*gammaOA          
        x30=xka/x
        x31= xka*xka2/x/x
        x38= x31*xka3/x
        
c        write(6,'(A,16E15.6)') 'aaaa', xkacid1,
c     &       xka,xka2, ml(29), ml(30), ml(31),gammah, gammahOA

        ML(29) = xmOA /(1 + x30 + x31+x38)
        ML(30) = ML(29)*x30         !oxalic acid ions
        ML(31)=ML(29) * x31 ! OA-2
        ML(38)=ML(29) * x38 ! A-3

        
        
        
        xka=10d0**(-3.84d0)/gammaH/gammaLA ! take the activity coefficient of H+
        aa=ml(33)+ml(34)
        ml(33) = aa/(1+ xka/x)
        ml(34)=  aa-ml(33)

        
           fcnNEW=0
        do j=1,20
           FCNNEW=fcnnew+ IZC(J)* ML(J)
           ENDDO

        do j=23,NP
           FCNNEW=fcnnew+ IZC(J)* ML(J)
           ENDDO
 
        fcnnew=-fcnnew
 

 22      continue

c     update gammah
         
         
      

      enddo
      
      
c      print*,'267', ml(27),ml(28)

      
 	return
	end

      
      
       subroutine vapnew(T0,M0,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poA,pl)
        IMPLICIT REAL*8 (A-H,O-Z)
        parameter (np=50,NSMM=100)         ! number of species
       real*8 gamma2(NSMM*2,NP)
      common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     &     , gammah2oA

      common /gamma2/gamma2
      common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

        real*8 m(NP)
        real*8 m0(NP)
        common /acid/ xkacid1, xkacid2,xhacid

      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
     & , imode_MA,imode_EDB
        
c        common/suls/ T,m
       common/suls/ T, M, ppartco2
      
       common /fvap/ alpha1,alpha0,ppvap,fvap_factor       
       real*8 x1(1),y1(1)
       
       data key_MA /0/
       real*8 ML_MA(100), gamma29_ma(100), gamma1_ma(100) 
       save key_MA, ML_MA, gamma29_ma, gamma1_ma
       character*100 path,filename
      common /path/path,Lpath
      logical ex
       common /ex/ ex
      
c     data key /0/
c       print*,fvap_factor
       
c     save  ff, key
        
c        if (imode_NH4NO3.eq.1) then
c        if (key.eq.0) then
c           key=1
c           ff=1
c           open(99,file='b.para')
c                      read(99,*,end=99) ff
c 99                   close(99)
c c       endif
c      else
c         ff=1
c      endif
      
           ff=fvap_factor

           
      do kk=1,np        
        M(kk)=M0(kk)
      enddo
      
        T=t0
       	call calHNew(T0,m0)
        do kk=1,np
        M(kk)=M0(kk)
        enddo
        call aw_back(t0,M0,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

        
        xh=gammaH
c     lactic acid
        xha = 2E6
        
       pL = 1013.5 * M0(33)/XhA
        
       xhh= 59* dexp( 4200*(1/T-1/298.15))
         pNH3 = m0(12)/m0(6)/ xknh3(T)*1013.5*gammaNH4/gammaH 
c         pNH3 = m0(23)*xhh*1013.5
         

         Xha =4000 * dexp(6200*(1/t-1/298.15d0)) ! Henrys law acetic acid M/bar
c     Sander 2015
         pacetic = 1013.5 * M0(9)/XhA

         xk= xkhx(T)*1000.d0**2/18.**2
c     Luo 1996
         phno3=m0(6)*m0(18)/xk*1013.5*gammaH*gammaNO3
c     Sanders 2015
         if (imode_NH4NO3.eq.1) then
            phno3=phno3/ff
            pnh3=pnh3/ff
         endif
         XHCO2 = 0.034 * dexp(2300 *(1/T-1/298.15))
  

         pco2 = m0(13) * 1013.5 /XHCO2*gamma2(Is,13)

         
         	xk=xhenry(T)

c     LUO 1996
         phcl=m0(6)*m0(17)/xk*gammaH*gammacl*1013d0

c     Compernolle and Müller 2014
          XhOA =6.1E8 * dexp(9800*(1/t-1/298.15d0)) ! Henrys law OA 
          if (xhacid.gt.0d0) xHOA = XHacid

c     ! oxalic acid

          poA = gammah2oA          *1013.5 * M0(29)/XhOA

          INQUIRE (FILE='AIOMFAC.dat', exist=ex)

          if (ex) then
             N11=11
             if (key_MA.eq.0) then
                key_MA=1

         filename='AIOMFAC.dat'
                open(99, file=filename)
                read(99,'(A1)')
                
                DO I=1,N11
                   read(99,*) xx1,xx2,xx3,xx4,xx5,xx6,xx7
                   
                   xx= xx5/(1-xx5)* m0(1)
                   print*,xx, xx7
                   mL_Ma(i)= XX
                   gamma29_ma(I)= xx7
                   write(55,*) xx, xx5*xx7

                enddo

                read(99,'(A1)')
                DO I=1,N11
                   read(99,*) xx1,xx2,xx3,xx4,xx5,xx6,xx7
c                  xx= xx5/(1-xx5)* m0(1)
c                   print*,xx, xx7
c                   mL_Ma(i)= XX
                   gamma1_ma(I)= xx7
                   write(55,*) xx, xx5*xx7

                enddo

c

             endif

             N1=1
             x1(1)= m0(29)
             
             call intpl(ml_MA, gamma29_MA,n11,x1,y1,n1)
             gamma29 = y1(1)
             call intpl(ml_MA, gamma1_MA,n11,x1,y1,n1)
             gamma1= y1(1)
             xmm=m0(1)+m0(2)

             DO I=6,np-npsolid
                xmm=xmm+m0(I)
             enddo

             pOA = XHOA * gamma29 * m0(29)/xmm ! poa
c             print*, xhOA, gamma29,m0(29), xmm
             
c             print*,'poa',poa
c     calculate the water activity of maleic acid

             
             aw=aw * m0(1)/(m0(1)+m0(29))* gamma1
          else
c     ideal for water
c             aw=aw * m0(1)/(m0(1)+m0(29))
          endif
          
C          stop
          
        return
        end
      


c     --------------------------------------------
c
c     Calculate the H+ and OH- concentration and dissociation  to ensure to charge balance 
c     
c     ------------------------------------------




c     -----------------------------------------------------------------------------
c
c     calculate the H+ and OH and discociations, CO2 is in Equilibrium with gas phase
c     defined by iseqCO2=1
c     
c     -----------------------------------------------------------------------------
      subroutine getm5(Ta,ma,ppartco2)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50,NSMM=100)
        real*8 fcnm5
	external  fcnm5
        real*8 MA(NP),m(NP)
	common/suls/ T, M, pp
        
       
        
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq



        common /kout/xx2

      common /output/imode_output,idiff,imode_pH,imode_eq

        real*8 xn(NSMM, np)
        common /xn/xn
       common /time/time, xn_area


       iseqco2=0
       pp=ppartco2
       t=ta
               do kk=1,np
        M(kk)= MA(kk)
      enddo
      
        call vapnew(T,M,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
         xmax = m(5)
         xmin=m(5)
      if (pco2.ge.pp) then
 17      xmin =xmin*.9 
         xx1=fcnm5(xmin)       ! pp -pco2vap
c         print*, 'xmin', xmin, xx1
         if (xx1.le. 0.001) goto 17
       
      endif
c      print*, xmin, xx1
      
            if (pco2.le.pp) then
 18      xmax =xmax*1.1 
         xx1=fcnm5(xmax)       ! pp -pco2vap
c         print*, 'xmax', xmax, xx1
         if (xx1.ge. 0.001) goto 18
         
      endif
         xx1=fcnm5(xmax)       ! pp -pco2vap
c      print*, xmax, xx1

	erabs=0.d0
	errel=0d0
	ITmax=200
         
      goto 333
      

        


 333    continue
        xx1=fcnm5(xmin)
        xx2=fcnm5(xmax)
c        print*,'xmin xx1 ', xmin,xx1
c        print*,'xmax xx2 ', xmax,xx2

	call dzbren3(fcnm5, erabs,errel,xmin,xmax,ITMAX)
 2      continue

        M(6)=xmax
        xx2= fcnm5(xmax)
        do kk=1,np
        MA(kk)=M(kk)
        enddo
        
 999            return
	end  

      subroutine calHNewco2(Ta,ma,ppartco2)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50,NSMM=100)
        real*8 fcnnew
	external  fcnnew
        real*8 MA(NP),m(NP)
	common/suls/ T, M,pp
        
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

        common /test/ Ierr



        common /kout/xx2
      common /output/imode_output,idiff,imode_pH,imode_eq


        real*8 xn(NSMM, np)
        common /xn/xn
       common /time/time, xn_area


       iseqco2=1
       pp=ppartco2


	t=ta
        do kk=1,np
        M(kk)= MA(kk)
        enddo

c        print*, 'T ' ,t

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	

c     H+ concentrations
        if (imode_output.eq.2 ) then
        mA(6)=1D-7
        mA(7)=1D-7
        return
        endif
        dd=.1D0

        if (is.le.0 .or. xn(Is,6).le.0d0) then
        xmin=1.E-14
 	xmax=10        ! 100% dissociation of acetic acid
        endif

        if (is.ge.1 .and. xn(Is,6).gt.0d0) then
           x6 = xn(is, 6)/xn(is,1)*M(1)
        
           xmin=x6/1.1
           xmax=x6*1.1
	xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)

        

        if (xx1*xx2.gt.0) then


        xx=10
        DO JJ=1,100
           xx=xx/2
	xx2= fcnnew(xx)
c        print*,xx,xx2
        if (xx2.ge.0d0) goto 223
        enddo
        print*,'error in calHCO2 !'
        stop
        
 223    continue
        xx=xx*2
        DO I=1,100
           xx=xx/1.1
	xx2= fcnnew(xx)
c        print*,xx,xx2
        if (xx2.ge.0d0) goto 224
      enddo
      
        print*,'error in calHCO2 !'
        stop
        
 224    xmin=xx
        xmax=xx* 1.1
c        print*, xx

        
              endif
           endif
           



        
	erabs=0.d0
	errel=0d0
	ITmax=200





	xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)






	if(xx1.le.0.) then
	xmax=xmin
	goto 2
	endif
	if(xx2.ge.0.) goto 2

	xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)
c           print*, 'xmin', xmin,xx1
c          print*, 'xmax', xmax,xx2
c        print*,'calhco2'
        
	call dzbren(fcnnew, erabs,errel,xmin,xmax,ITMAX)
 2      continue
c        print*, xmax
        M(6)=xmax
        xx2= fcnnew(xmax)
        do kk=1,np
        MA(kk)=M(kk)
      enddo
      
        if (ierr.eq.1) then
           print*, xmax, xx2
           print*,'error: calhco2'
           
           stop
        endif
        


        
 999            return
	end  
      


c     -------------------------------------------------------------
c
c     Calculate inactivaton times for SARS-CoV-2
c     
c     -------------------------------------------------------------


            subroutine cal_tau_sars(T,RH,pH,tau)
      implicit real*8 (a-h,o-z)
      parameter (N=3)
      real*8 pharr(N), tauarr(N), k(N),x(4),y(1),klog(N)
      real*8 x1(1),y1(1), ph5(111), tau5(111)
      data key /0/
      save key,x , ph5, tau5,n5
      


      if(key.eq.0) then
         ph5(1)=5.624
         ph5(2)=7.5
         ph5(3)=9.5
         ph5(4)=10
         ph5(5)=11
         ph5(6)=12
         ph5(7)=13
         
         tau5(1)=400950/dlog(100d0)
         tau5(2)=1D5/dlog(100d0)
         tau5(3)= 3700
         tau5(4)= 1150
         tau5(5)= 3190/dlog(100d0)
         tau5(6)= 117*2/dlog(100d0)
         tau5(7)= 20/dlog(100d0)
         
         N5=7
         do I=1,N5
            tau5(i)= dlog(tau5(I))
         enddo
         
      endif
      
      
      x(1)=0.50000E+01 
        x(2)=0.54638E+01 
        x(3)=0.24943E+01 
        x(4)=0.21971E+01 

          x(1)=      0.57834E+01
          x(2)=    0.50000E+01   
         x(3)=   0.20308E+01 
         x(4)=   0.21729E+01 

        
        phh=ph-x(4) 
      tau_sars = x(1)   + x(2)*datan(x(3)   *phh )
      tau_sars =  dexp(tau_sars)
      tau = tau_sars/dlog(100d0)
c     base

      if (ph.ge.ph5(1)) then
            x1(1)=ph

            N1=1


            
           call intpl(ph5, tau5, n5, x1,y1, n1)
           
          tau= dexp(y1(1))
c          tau = tau5(1) + (tau5(4)-tau5(1))/( ph5(4)-ph5(1))*(ph-ph5(1))
c          tau= dexp(tau)
          
c          print*, ph,tau
          
          
       endif
         
      
      
      return
      end




c     -------------------------------------------------------------
c
c     ions diffusion coefficent, obtained from EBD data, scaled to 293 K 
c     
c     -------------------------------------------------------------

      subroutine cal_ions(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)
      
      common /enhance/radius,r1,enh_factor
c      common /Ienhance/Ienh,iscenter
      common /Ienhance/Ienh,iscenter,isahoo
      data key /0/
      data fac1 /1d0/
      data fac2 /1d0/
      logical ex
      save fac1,fac2 , key,ex,fac3,fac4
      if (key.eq.0) then
      INQUIRE (FILE='DL_SLF.dat', exist=ex)
      fac1=1d0
      fac2=1d0
      fac3=1d0
      fac4=1d0
      key=1
      
      if ( ex ) then
         open(99, file='DL_SLF.dat')
         read(99,*) fac1
         read(99,*) fac2
         read(99,*,end=23) fac3
         read(99,*,end=23) fac4

c         print*, fac1,fac2
c         stop
         

 23      continue

         close(99)
      endif
      endif
      

      aw1=1

      t0=298.15
      call cal_dlaw_citric(T0,aw1,dl298)
      t293=293.15
      call cal_dlaw_citric(T293,aw1,dl293)
c      print*, dl0
      
      Dl0= dl293/dl298 * 2.44E-5  ! scale d0 to 2.44E-5 cm2/s at T0 and aw=1

      enhance = 1


      if (r1.ge.1D-12 .and. Ienh.eq.1) then

         xlc = 340d-4/20d0**3*(radius*1D4)**3

         enhance = xlc/2/r1
c         print*, enhance, radius, xlc,r1

         
      if (enhance.lt.1d0) enhance=1d0
         endif

      enhance = 1

         ff=0.54*(enhance)
            ff1=1/1.46

            fac(1)=.06D-9*ff*ff1*fac1
            fac(2)=0.07E-9*ff*ff1*fac2
            fac(3)=0.15E-9*ff*ff1*fac3
            fac(4)= .4D-9*ff *.9*ff1*fac4
            fac(5)= 1.3D-9*ff*1.1*ff1
       fac(6)= 1.8E-8*ff**(1/2d0)*ff1**(0.5d0)
       fac(7)= 6E-8 
       fac(8)= dl0/10
       fac(9)= dl0/2
       fac(10)= dl0
       
      rh(1)= .0
      rh(2)=.5
      rh(3)=.575
      rh(4)=.65
      rh(5)=.72
      
       rh(6)=0.85
       rh(7)=0.90
       rh(8)=0.96
       rh(9)=0.99
       rh(10)=1
       DO kk=1,10
      fac(kk)=dlog(fac(kk))
      enddo
      N1=1
       N5=10
       
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))

c     no reduction of diffusion coefficient for crystal in the center
       if (iscenter.eq.1) return
c     DO reductions for co-shells crystals
       aww=aw0
       if( aww.le.0.5d0) Aww=0.5
       ff=dlog(0.4d0)*(1-aww)/0.5d0
       dl = dexp(ff)* dl        ! multiply by 0.4 at aw <= 0.5
       
       
       return
       end



c     --------------------------------------------------
c
c    activity product crystal = a_Na * a_Cl
c
c     --------------------------------------------------

      
      function apNaCL(t)
      IMPLICIT REAL*8 (A-H,O-Z)
c          apnacl=      0.13852E+02 -0.21891D-1 *(t-298.15)
c      apnacl= dexp(0.30606E1+ 129.22/t) !pinho 2005
c     apnacl= dexp(3.2635+ 132.86/t) !pinho 2005
       apnacl= dexp(3.4734+ 166.89/t) !pinho 2005      
c       apnacl= dexp(4.5473+ 188.3/t) !pinho 2005      

      return
      end
      


c     --------------------------------------------------
c
c     calculate the flux f1 (from the shell to to crystal) and f2 (the flux liquid to crystal of the same shell)
c
c     NS: the total number of shells
c      II: the location of the crystal 
c     
c     --------------------------------------------------

      subroutine  cal_Misch_flux (NS, iI,Tdrop, x,f1,f2)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      parameter (np=50,NSMM=100,npsolid=12)
      real*8  x(*),ml(NP),xa(NSMM),ml0(np),ml1(NSMM)
      
       common /eutectic/ feutectic0
      real*8 dl_factor2(2,np),dl_factor(NP)
c     
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin

      real*8 xn(NSMM,NP),x1(NP)
      common /xn/xn
      


      common /output/imode_output,idiff,imode_pH,imode_eq
      
     
      common /Ienhance/Ienh,iscenter,isahoo
      common /pi/pi
      


      real*8 MM(NP) ! molar Mass
      real*8 Mv(NP) ! molar volume
      common /M/ MM,mv          ,izc(np)
      common /Nacl/T,ML,ML0
      external fcn_ap

      common /testA/Ierr
       feutectic=feutectic0


       
      Ierr=0

      I=II

      f1=0d0

      xvsolid=mv(Np) *xn(I,np) + mv(Np-1) *xn(I,np-1)
      
      vol= 4*pi/3d0*x(I)**3 + xvsolid
      xvliquid =4*pi/3d0*( x(I+1)**3-x(I)**3) -xvsolid

      r1=(vol/4d0/pi*3)**(1/3d0)
      r2=x(I+1)

      if (I.lt.NS) then
         if (r2.le. x(I+2)/2) r2 = x(I+2)/2
      endif
      
      v1= xvliquid
      
            
      ML(1) = 1000d0/mm(1)

      t=Tdrop

c     take the mean composition of shell I and I+1
      DO J=1,NP
         ML1(J)=ML(1)* xn(I,J)/xn(I,1)
                  if (I.lt.NS) then
         ML(j)= (ML(1)* xn(I,J)/xn(I,1)+ ML(1)* xn(I+1,J)/xn(I+1,1))/2
         
        else
         ML(j)= ML(1)* xn(I,J)/xn(I,1)
         endif
c         print*, i,j, ML(j) 
      enddo
        do kk=1,np
           ML0(kk)=ML(kk)
           enddo
      xmax=50d0

       xmin=1D-3
       xmax =xmax*.999999


       erabs=0.
	errel=0d0
        itmax=100
c        pi=dacos(-1d0)

      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
      Ienh=0
      if (iscenter.eq.1) ienh=1
      call cal_dlaw(T,aw,dl)
c      if (idiff.eq.3)   call cal_dlaw_suc(T,aw,dl)
      if (idiff.eq.4)   call cal_dlaw_citric(T,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
      if (idiff.eq.6)   call cal_dlaw_walker_h2o(T,aw,dl)

       dl_factor(16)=dl_factor2(1,16)
       dl_factor(17)=dl_factor2(1,17)
       dl_factor(11)=dl_factor2(1,11)
       dlfNa=dl_factor(16)
       dlfCl=dl_factor(17)

c       print*, 'before ap ', ml0(16),ml0(17)
       

c        print*, 'cal_mish befor', ierr
       xmin=1D-3
       xmax =xmax*.999999
c       ymin=fcn_ap
c       ymax=fcn_ap
       
c       print*,'misch_flux'
       call dzbren(fcn_ap,erabs,errel,xmin,xmax,ITMAX)
        

        if (Ierr.ge.1) then
           print*, 'erorr: cal_misch 2'
           stop
           endif

c       print*, 'after ap ', ml0(16),ml0(17)
c        call calHNew(T,mL0)
       xx2=fcn_ap(xmax)

       dlna= dl* dl_factor(16)
       DlCl= dl* dl_factor(17)
       Dl11= dl* dl_factor(11)

       c1 = xn(I,1)/v1* ML(16)/ML(1)
       c0 = xn(I,1)/v1* ML0(16)/ML0(1)
       
       fNa = dlna  * (c1-c0)*4 *pi*r1*(r2)/(r2-r1)
       fna=fna 
       
       c1 = xn(I,1)/v1* ML(17)/ML(1)
       c0 = xn(I,1)/v1* ML0(17)/ML0(1)
       fcl= dlcl  * (c1-c0)*4*pi*r1*(r2)/(r2-r1)
       fcl=fcl 


c       print*, fcl, fna
       ff = 1 !vliquid/VOL
       if (I.eq.1) ff=1d0
       f2=fcl*ff
c 
       
       call aw_back(t,ML1,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
       Amisch1= ML1(16)*ml1(17)*gammaCl*gammaNa 

       if (Amisch1.le.apnacl(T) .and. f2 .gt.0d0) f2=0d0 !  no growth when the I shell is not saturated
       if (Amisch1.ge.apnacl(T) .and. f2 .lt.0d0) f2=0d0 !  no dissolve when NaCl is super-saturated 
       
       
       if (I.ge.2) then
          
c     calculate f1
c     the flux from inner shell to crystal  !!!!!
c     c0-c1
      r1=x(I)  ! outer shell liquid I-1
      r2=x(I-1) ! inner radius of shell I-1 liquid

       v2 = 
     & 4*pi/3 * x(I-1)**3 + mv(Np)*xn(I-1,np) + mv(Np-1)*xn(I-1,np-1)
      r2=dsqrt(v2/4/pi*3)  ! inner radius of shell I-1


      vliq= xn(I-1,1) *mv(1)+ xn(I-1,2)*mv(2)
      DO J= 6,NP-npsolid
         vliq= vliq + xn(I-1,J)*mv(J)
      enddo
      
         
c       ff =vsolid/(vliq+vsolid)
c     take the mean compostion of the shell I and I+1
      DO J=1,NP
         I1= I-1
         ML1(j)= ML(1)* xn(I1,J)/xn(I1,1)
         ML(j)= (ML(1)* xn(I,J)/xn(I,1)+ ML(1)* xn(I1,J)/xn(I1,1))/2
      enddo

        do kk=1,np
      
      ML0(kk)=ML(kk)
      enddo
      



      xmax=100d0

      xmin=1D-3
       xmax =xmax*.999999

       xx1=fcn_ap(xmin)



       erabs=0.
	errel=0d0
        itmax=100
c        pi=dacos(-1d0)
      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        


            Ienh=0
      if (iscenter.eq.1) ienh=1

      call cal_dlaw(T,aw,dl)

c      if (idiff.eq.3)   call cal_dlaw_suc(T,aw,dl)
      if (idiff.eq.4)   call cal_dlaw_citric(T,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
      if (idiff.eq.6)   call cal_dlaw_walker_h2o(T,aw,dl)

       dlfNa=dl_factor2(1,16)
       dlfCl=dl_factor2(1,17)

c       print*,'mish_flux1'
       call dzbrens(fcn_ap,erabs,errel,xmin,xmax,ITMAX)

        if (ierr.ge.1) then
           print*, 'error: cal_misch 1'

           stop
           endif

       xx2=fcn_ap(xmax)

       dlna= dl* dl_factor2(1,16)
       DlCl= dl* dl_factor2(1,17)

       c1 = xn(1,1)/vliq* ML(16)/ML(1)
       c0 = xn(1,1)/vliq* ML0(16)/ML0(1)
       
       fNa = dlna  * (c1-c0)*4 *pi*r1*(r2)/(r1-r2)

       c1 = xn(1,1)/vliq* ML(17)/ML(1)
       c0 = xn(1,1)/vliq* ML0(17)/ML0(1)
       fcl= dlcl  * (c1-c0)*4*pi*r1*(r2)/(r1-r2)
c     oppsite sign as F2

       f1=(fcl+fna)/2 ! *vsolid/(v1+vsolid)/2


       call aw_back(t,ML1,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
       Amisch1= ML1(16)*ml1(17)*gammaCl*gammaNa 
c     f1 and f2 has oppsote sign
       
       if (Amisch1.le.apnacl(T) .and. f1 .gt.0d0) f1=0d0 !    nogrowth when the I shell is not saturated
       if (Amisch1.ge.apnacl(T) .and. f1 .lt.0d0) f1=0d0 !  no desolve when super saturated

          

       endif

c     introduce enhance memnta factor
c
c       thinmin = 1D-7 ! 

       
       return
        end
      


c     --------------------------------------------------
c
c     Required by cal_misch_flux, to calculate the equilibrium Na and Cl concentration over NaCl Crystal
c      
c     --------------------------------------------------

      function fcn_ap(xmcl)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      parameter (np=50)
      real*8 ml(NP),ml0(NP)

     

       
      real*8 dl_factor2(2,np),dl_factor(NP)
     
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin

      common /Nacl/T,ML,ML0


      common /output/imode_output,idiff,imode_pH,imode_eq

  


      
c      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        


       ML0(17)=xmcl
       d17=dl_factor2(1,17)
       if (d17 .le.1D-30) d17=1d0
       d16=dl_factor2(1,16)
       if (d16 .le.1D-30) d16=1d0

       ML0(16)= d17/d16 * (ML0(17)-ML(17))+ ML(16)
      
       if (ML0(16).le.1D-10 )  ml0(16)=1D-10


       call aw_back(t,ML0,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
 
c       print*, 'ml0 ',ml0(16),ml0(17) 
       Amisch1= ML0(16)*ml0(17)*gammaCl*gammaNa 
       fcn_ap= Amisch1-Apnacl(T)
       
       return
      end
      


      function fcnm5(x)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      parameter (np=50)
      real*8 m(NP),m0(NP),fcnm5
	common/suls/ T, M,pp
      
c      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        

        m(5)=x
        
      call vapnew(T,M,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
 
c       print*, 'ml0 ',ml0(16),ml0(17) 

      fcnm5= (pp-pco2)/pp

      

       
       return
      end
      

      
c     ------------------------------------------------
c
c     calculate H2O diffusion coefficient in Pinene solution Lienhard et al
c
c     ------------------------------------------------

      subroutine       cal_dlaw_piene(T,aw,dl,xw)
      implicit real*8 (a-z)
      aw1=1d0
      call  cal_dlaw_citric(t,aw1,dw1)
      dw0 = 7E-11* dexp ( -65500/8.314*(1/T-1/300d0))
      alpha=1
      ta=T
      if (T.ge.273d0) ta =273d0
      A= -18.31+0.063*TA
      B= -10.65 + .039*TA
      alpha= (1-xw)**2 *( A +3*b-4*b*(1-xw))
c      print*,xw,dw1,alpha
      alpha=dexp(alpha)
      
      Dl = dw0** (1-xw*alpha) * dw1**(xw*alpha)
      return
      end
      

c     -------------------------------------------------------
c
c     calculate aw Pinene solution Lienhard et al
c
c     -------------------------------------------------------
      
      function       aw_pinene(xw)
      implicit real*8 (a-z)
      xmfs  = (1-xw)* 150 / (xw*18 + (1-xw)* 150)
      aw_pinene = (1-xmfs)/( 1 - 0.85848*xmfs-0.09026*xmfs**2)

      return
      end

c     ---------------------------------
c
c     calculate the OH- activity coefficient in a NaOH solution
c
c     ---------------------------------
      
      function get_gammaOH(xMplus)
      implicit real*8 (a-h,o-z)
      parameter (N=240,n1=1)
      real*8  x(N),y(N),x1(N1),y1(N1)
      data key /0/
      data x /
     &    0.0000,     0.0010,     0.0020,     0.0030,     0.0040, 
     &    0.0050,     0.0060,     0.0070,     0.0080,     0.0090, 
     &    0.0100,     0.0110,     0.0120,     0.0130,     0.0140, 
     &    0.0150,     0.0160,     0.0170,     0.0180,     0.0190, 
     &    0.0200,     0.0210,     0.0220,     0.0230,     0.0240, 
     &    0.0250,     0.0260,     0.0270,     0.0280,     0.0290, 
     &    0.0300,     0.0310,     0.0320,     0.0330,     0.0340, 
     &    0.0350,     0.0360,     0.0370,     0.0380,     0.0390, 
     &    0.0400,     0.0410,     0.0420,     0.0430,     0.0440, 
     &    0.0450,     0.0460,     0.0470,     0.0480,     0.0490, 
     &    0.0500,     0.0510,     0.0520,     0.0530,     0.0540, 
     &    0.0550,     0.0560,     0.0570,     0.0580,     0.0590, 
     &    0.0600,     0.0610,     0.0620,     0.0630,     0.0640, 
     &    0.0650,     0.0660,     0.0670,     0.0680,     0.0690, 
     &    0.0700,     0.0710,     0.0720,     0.0730,     0.0740, 
     &    0.0750,     0.0760,     0.0770,     0.0780,     0.0790, 
     &    0.0800,     0.0810,     0.0820,     0.0830,     0.0840, 
     &    0.0850,     0.0860,     0.0870,     0.0880,     0.0890, 
     &    0.0900,     0.0910,     0.0920,     0.0930,     0.0940, 
     &    0.0950,     0.0960,     0.0970,     0.0980,     0.0990, 
     &    0.1000,     0.1100,     0.1200,     0.1300,     0.1400, 
     &    0.1500,     0.1600,     0.1700,     0.1800,     0.1900, 
     &    0.2000,     0.2100,     0.2200,     0.2300,     0.2400, 
     &    0.2500,     0.2600,     0.2700,     0.2800,     0.2900, 
     &    0.3000,     0.3100,     0.3200,     0.3300,     0.3400, 
     &    0.3500,     0.3600,     0.3700,     0.3800,     0.3900, 
     &    0.4000,     0.4100,     0.4200,     0.4300,     0.4400, 
     &    0.4500,     0.4600,     0.4700,     0.4800,     0.4900, 
     &    0.5000,     0.5100,     0.5200,     0.5300,     0.5400, 
     &    0.5500,     0.5600,     0.5700,     0.5800,     0.5900, 
     &    0.6000,     0.6100,     0.6200,     0.6300,     0.6400, 
     &    0.6500,     0.6600,     0.6700,     0.6800,     0.6900, 
     &    0.7000,     0.7100,     0.7200,     0.7300,     0.7400, 
     &    0.7500,     0.7600,     0.7700,     0.7800,     0.7900, 
     &    0.8000,     0.8100,     0.8200,     0.8300,     0.8400, 
     &    0.8500,     0.8600,     0.8700,     0.8800,     0.8900, 
     &    0.9000,     0.9100,     0.9200,     0.9300,     0.9400, 
     &    0.9500,     0.9600,     0.9700,     0.9800,     0.9900, 
     &    1.0000,     1.2000,     2.2000,     3.2000,     4.2000, 
     &    5.2000,     6.2000,     7.2000,     8.2000,     9.2000, 
     &   10.2000,    11.2000,    12.2000,    13.2000,    14.2000, 
     &   15.2000,    16.2000,    17.2000,    18.2000,    19.2000, 
     &   20.2000,    21.2000,    22.2000,    23.2000,    24.2000, 
     &   25.2000,    26.2000,    27.2000,    28.2000,    29.2000, 
     &   30.2000,    31.2000,    32.2000,    33.2000,    34.2000, 
     &   35.2000,    36.2000,    37.2000,    38.2000,    39.2000, 
     &   40.2000,    41.2000,    42.2000,    43.2000,    44.2000, 
     &   45.2000,    46.2000,    47.2000,    48.2000,    49.2000/
      data y/
     &    1.0000,     0.9593,     0.9407,     0.9259,     0.9131, 
     &    0.9016,     0.8911,     0.8813,     0.8721,     0.8635, 
     &    0.8553,     0.8474,     0.8399,     0.8327,     0.8257, 
     &    0.8190,     0.8125,     0.8061,     0.8000,     0.7941, 
     &    0.7883,     0.7826,     0.7771,     0.7718,     0.7665, 
     &    0.7614,     0.7564,     0.7515,     0.7467,     0.7420, 
     &    0.7374,     0.7328,     0.7284,     0.7240,     0.7198, 
     &    0.7156,     0.7114,     0.7074,     0.7034,     0.6995, 
     &    0.6956,     0.6918,     0.6881,     0.6844,     0.6807, 
     &    0.6772,     0.6737,     0.6702,     0.6668,     0.6634, 
     &    0.6601,     0.6568,     0.6535,     0.6503,     0.6472, 
     &    0.6441,     0.6410,     0.6380,     0.6350,     0.6320, 
     &    0.6291,     0.6262,     0.6234,     0.6206,     0.6178, 
     &    0.6150,     0.6123,     0.6096,     0.6070,     0.6043, 
     &    0.6017,     0.5992,     0.5966,     0.5941,     0.5916, 
     &    0.5892,     0.5867,     0.5843,     0.5819,     0.5796, 
     &    0.5772,     0.5749,     0.5726,     0.5704,     0.5681, 
     &    0.5659,     0.5637,     0.5615,     0.5594,     0.5573, 
     &    0.5551,     0.5531,     0.5510,     0.5489,     0.5469, 
     &    0.5449,     0.5429,     0.5409,     0.5389,     0.5370, 
     &    0.5351,     0.5167,     0.4998,     0.4842,     0.4698, 
     &    0.4563,     0.4438,     0.4321,     0.4212,     0.4108, 
     &    0.4011,     0.3920,     0.3834,     0.3752,     0.3674, 
     &    0.3601,     0.3531,     0.3465,     0.3402,     0.3341, 
     &    0.3284,     0.3229,     0.3176,     0.3126,     0.3078, 
     &    0.3032,     0.2988,     0.2945,     0.2904,     0.2865, 
     &    0.2827,     0.2791,     0.2756,     0.2722,     0.2690, 
     &    0.2659,     0.2628,     0.2599,     0.2571,     0.2544, 
     &    0.2518,     0.2492,     0.2468,     0.2444,     0.2421, 
     &    0.2399,     0.2377,     0.2356,     0.2336,     0.2316, 
     &    0.2297,     0.2279,     0.2261,     0.2243,     0.2227, 
     &    0.2210,     0.2194,     0.2179,     0.2164,     0.2149, 
     &    0.2135,     0.2121,     0.2108,     0.2095,     0.2082, 
     &    0.2069,     0.2057,     0.2046,     0.2034,     0.2023, 
     &    0.2013,     0.2002,     0.1992,     0.1982,     0.1972, 
     &    0.1963,     0.1954,     0.1945,     0.1936,     0.1927, 
     &    0.1919,     0.1911,     0.1903,     0.1895,     0.1888, 
     &    0.1881,     0.1874,     0.1867,     0.1860,     0.1853, 
     &    0.1847,     0.1749,     0.1684,     0.1909,     0.2288, 
     &    0.2803,     0.3468,     0.4305,     0.5349,     0.6639, 
     &    0.8224,     1.0159,     1.2508,     1.5341,     1.8740, 
     &    2.2790,     2.7586,     3.3226,     3.9815,     4.7456, 
     &    5.6254,     6.6308,     7.7710,     9.0538,    10.4853, 
     &   12.0693,    13.8072,    15.6969,    17.7328,    19.9053, 
     &   22.2005,    24.6002,    27.0816,    29.6178,    32.1779, 
     &   34.7276,    37.2297,    39.6452,    41.9340,    44.0561, 
     &   45.9729,    47.6480,    49.0486,    50.1465,    50.9192, 
     &   51.3502,    51.4303,    51.1570,    50.5356,    49.5783/

      save x ,y,key
      character*100 path,filename
      common /path/path,Lpath

      if (key.eq.0) then
         key=1
         filename=path(1:lpath)//'IVEA/model/NaOH.g'
         open(99, file=filename)
         DO I=1,N
            read(99,*) x(I),y(I)
         enddo
         close(99)
      endif
      
         

      

      x1(1)= xmplus
      call intpl(x,y,n,x1,y1,n1)
      get_gammaOH=y1(1)
      return
      end

      
c     --------------------------------------------
c      
c     calcalate inactivation time of IAV virus is second
c
c     --------------------------------------------

      function tau_ivea(ph,xm0)
      implicit real*8 (a-h,o-z)
      
      real*8 ph, tau_ivea, tau_ivea_1,xx,dd,xm0,xm,x(16), phh
      integer i,ii,key
      real*8 phb(10),taub(10),fac75(10),x1(1),y1(1)
      
      
      data key /0/
      save key,x,phb,fac75,nb
       

      if(key.eq.0) then
         key=1
        x(1) =   0.50131E+01 
           x(2)=    0.62214E+01   
          x( 3)=         0.12972E+01   
          x(4)=    0.53802E+01   
          x(5)=    0.14134E+01   
          x(6)=   -0.47140E+00   
          x(7)=   -0.20450E+01   
          x(8)=    0.46916E+00   

          NB=4
         ph7=7.5
          phb(1)=ph7
          phb(2)=9
          phb(3)=10
          phb(4)=12

          taub(2)=1.1E6/dlog(100d0)
          taub(3)=1.E4/dlog(100d0)
          taub(4)=1.6E2/dlog(100d0)


         xm00=0d0
         tau75=tau_ivea1(ph7,xm00)
         fac75(1)=0d0
         taub(1)=tau75
         

          DO I=1,NB

        
         fac75(I) = dlog(taub(I)/tau75)
          print*,'phh' , phb(I), fac75(I)
c         print*,phh, tau_IVEA, tau75
         
          enddo
          

       close(25)
      endif
      
       xm=xm0-0.027
       if (xm.le.0d0) xm=0d0

         phh=ph-x(4) + xm*x(5)
         if (phh.ge.7.5d0) phh= 7.5
         tau_ivea = x(1)   + x(2)*datan(x(3)*(1+xm*x(6))   *phh )

         taul = x(7)+ x(8)* ph
         tau_ivea = ( dexp(tau_ivea) + dexp(taul))

         
          if (ph.ge.7.5) then
         n1=1
         x1(1)=ph
         call intpl(phb, fac75,NB,x1,y1,n1)

         tau_ivea =tau_ivea *dexp(y1(n1))

         endif


             return
      end 
      

      
      function tau_ivea1(ph,xm0)
      implicit real*8 (a-h,o-z)
      
      real*8 ph, tau_ivea, tau_ivea_1,xx,dd,xm0,xm,x(16), phh
      integer i,ii,key
      real*8 phb(10),taub(10),fac75(10),x1(1),y1(1)
      
      
      data key /0/
      save key,x,phb,fac75,nb
       

      if(key.eq.0) then
         key=1
        x(1) =   0.50131E+01 
           x(2)=    0.62214E+01   
          x( 3)=         0.12972E+01   
          x(4)=    0.53802E+01   
          x(5)=    0.14134E+01   
          x(6)=   -0.47140E+00   
          x(7)=   -0.20450E+01   
          x(8)=    0.46916E+00   

          NB=4
          phb(1)=7
          phb(2)=9
          phb(3)=10
          phb(4)=12

          taub(2)=1.5E6/dlog(100d0)
          taub(3)=1.E4/dlog(100d0)
          taub(4)=1.6E2/dlog(100d0)

             phh=7.5
         tau_ivea = x(1)   + x(2)*datan(x(3)   *phh )
         taul = x(7)+ x(8)* phh
         tau_ivea = ( dexp(tau_ivea) + dexp(taul))
         tau75=tau_ivea
         fac75(1)=0d0

          DO I=2,NB
             phh=phb(I)
        
         fac75(I) = dlog(taub(I)/tau75)
c         print*,'phh' , phh, fac75(I)
c         print*,phh, tau_IVEA, tau75
         
          enddo
          

       close(25)
      endif
      
       xm=xm0-0.027
       if (xm.le.0d0) xm=0d0

         phh=ph-x(4) + xm*x(5)
         if (phh.ge.7.5d0) phh= 7.5
         tau_ivea = x(1)   + x(2)*datan(x(3)*(1+xm*x(6))   *phh )

         taul = x(7)+ x(8)* ph
         tau_ivea1 = ( dexp(tau_ivea) + dexp(taul))

         
          if (ph.ge.17.5) then
         n1=1
         x1(1)=ph
         call intpl(phb, fac75,NB,x1,x1,n1)
         tau_ivea =tau_ivea *dexp(y1(1))

         endif


             return
      end 

      

c     ---------------------------------------------------------
c
c     calculate the molality of the solution at given aw (in equlibrium) and T
c
c     ---------------------------------------------------------

      subroutine cal_ml(T,rh,ML)

         
	 implicit real*8 (a-h, o-z)
	 external  fml

	parameter (NP=50)
	real*8 mm(NP),mv(NP),mL(NP),mlA(NP),ml0(NP) ! mole mass
	common /m/ mm, mv,izc(np)
	common /mLA/ rh1,t1, ml0,MLA
        common /testA/Ierr

        Ierr=0

        do kk=1,np
         
           ML0(kk)=ML(kk)
           enddo
c         pwrint*, 'x', x, MV(1)

         xms=MM(NP)*ML(NP)+MM(2)*ML(2)
         DO I=6, np
            xms=xms+ MM(I)*ML(I)
            enddo
            xmfs0=xms/1000


	 t1=t
	 rh1=rh

c     x = solute /H2O ratio
	 xmin=.5
	 xmax=1.5
         

	 erabs=0.d0
	 errel=0d0
	 ITmax=1000

c         xmin = 1D-10
 	 ymin=fml(xmin)
	 ymax=fml(xmax)
c	 ymax1=fml(15d0)
 
c         print*,'xmin', xmin,ymin
c         print*,'xmax', xmax,ymax
         
         if (ymin.ge.0) then

 	 ymin=fml(xmin)
            goto 12
            endif
            if (ymax.le.0) goto 12

c         print *, 'xml' , xmin,xmax
c         print *, 'fml' , ymin,ymax

	 call dzbrens(fml,erabs,errel,xmin,xmax,ITMAX)


         if (Ierr.ge.1) then
            print*, 'error: Cal_ml'
            stop
            endif
	 ymax=fml(xmax)

  12	 continue
         do kk=1,np
            ML(kk)=MLA(kk)
            enddo
c         call calhnew(t,mla)
      call aw_back
     & (t,MLA,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
 
  
	 return

	 end

c     ---------------------------------------------------------
c
c     required ny cal_ml, to calculate the molalities at given aw
c      
c     ---------------------------------------------------------
      
	 function fml(X) 
	 implicit real*8 (a-h,o-z)
	parameter (NP=50)
	 real*8 mm(NP),mv(NP),mL0(NP),mla(NP) ! mole mass

	common /m/ mm, mv,izc(NP)
	common /mLA/ rh,t, ml0,MLA

c        print*, x,xmfs0



          MLA(1)= 1000d0/mm(1)
           DO I=2,NP
            MLA(I)= ML0(I)*x
           enddo
           
      call calhnew (t,mlA)
      call aw_back
     & (t,MLA,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

 	   Fml=(RH-aw)

c           print*, 'fml', x, aw
           


 10 	return
	 end
 
c     ---------------------------------------------------------------------------------------
c      
c     calculate water activity and acitivity coeffcients of ions using Piter ion interaction model for ions
c      
c     for neutral species, The Rauolts Law is used
c
c     ---------------------------------------------------------------------------------------
 
      subroutine aw_back_model
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
        implicit real*8 (a-h,m,o-z)
      integer NP
      parameter (np=50,npsolid=12,NSMM=100)
      real*8 ML(Np)

      real*8 MM(NP) ! molar mass
      real*8 Mv(NP) ! molar volume
      integer IZC(np)
      common /M/ MM,mv         ,izc 
c      common /awin/ awin     ,aws   ,xvol
      common /awin/ awin,aws,xvol,xmi
      
	Parameter( NMAx=12 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax)
 	real*8 C1(nmax,nmax),omega(nmax,nmax),mc0(nmax),ma0(nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax),xs(100),xfit(10)

 	real*8 b0z(nmax,nmax),B1z(nmax,nmax),C0z(nmax,nmax)
 	real*8 C1z(nmax,nmax),omegaz(nmax,nmax)

        
	integer Iflag(NMAX,NMAX)
	data ZC /12*1./,ZA/12*1./
        common /a/ ah2so4



      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
     & , imode_MA,imode_EDB

      common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     &     , gammah2oA,gammaPO4
        

        DATA KEYA /0/
	real*8 B20(nmax,nmax),alpha2(nmax,nmax)
	common /alpha/ alpha2, b20

c                      read(99,*) ii,
c              read(99,*) ii,
c              read(99,*) ii,
c              read(99,*) ii,xlamnn
      parameter (ndim=35)

      
      
      
       
      common/lam/xlamNa,xlamcl
        common /xm24/ xm24,xlamc,xlama1,xlama2,xlamnn
        common /xm29/ xm29,xlam1_29,xme1_29,xlam2_29,xme2_29

              

        
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

        real*8 gamma2(NSMM*2,NP)
        common /gamma2/gamma2

        
        data xfit /
     1    0.65165E-02, 
     2    0.92664E-01, 
     3    0.19114E-04, 
     4    0.78199E-01 ,
     5    0.12843E+01 ,
     6    0.79342E-02, 
     7   -0.72711E-01, 
     8    0.69578E-04 ,
     9    0.22108E-01 ,
     1    0.50531E+00 /
        save b0,B1,C0,c1,omega, keya,za,zc,
     &       b0z,B1z,C0z,c1z,omegaz,iflag
      character*100 path,filename
      common /path/path,Lpath

        
        if (is.le.1) is=1
        
        

        if (keya.eq.0) then
           keya=1
           DO I=1,nmax
           DO j=1,nmax
           omega(i,j)=1
           b0(i,j)=0d0
           b1(i,j)=0d0
           c0(i,j)=0d0
           c1(i,j)=0d0
           alpha2(i,j)=1d0
           b20(i,j)=0
           enddo
           enddo
           NC=4
           close(99)

c NaHCO3 ! Pitzer 1980
c     MA 1: HNO3
c     MA 2: HCl
c     MA 3: HSO4-
c     MA 4: SO42-
c     MA 5: CO3
c     MA 6: HCO3
c     MA 7: H2PO4-
c     MA 8: HPO42-
c     MA 9: HC2O4-
c     MA 10: C2O42-

           
c           NC
c     1: H+
c     2: NH4+
c     3: Na+
c     4: Mg++
c     5:K
c     6:Ca
            filename=path(1:Lpath)//'/IVEA/model/Pitzer.para'           
           open(99,file=filename)
        NC=6
        NA=10

           DO I=1,Nc
           DO j=1,NA
              read(99,*) ii,jj, b0(i,j), b1(i,j), c0(i,j),c1(i,j),
     & omega(i,j),b20(i,j),alpha2(i,j) 
        enddo
      enddo
c     replace Na-HPO4
          filename=path(1:Lpath)//'/IVEA/model/b_Na2HPO4.dat'           
           open(99,file=filename)

          I=3
          J=8
c          DO I =1,5
             read(99,*) ii ,b0(i,j)
             read(99,*) ii ,b1(i,j)
             read(99,*) ii ,c0(i,j)
             read(99,*) ii ,c1(i,j)
             read(99,*) ii ,omega(i,j)

      
          filename=path(1:Lpath)//'/IVEA/model/b_H3PO4_new'           
           open(99,file=filename)
        II= 1
        JJ= 7
              read(99,*) i, b0(II,JJ)
              read(99,*) i, b1(II,JJ)
              read(99,*) i, c0(II,JJ)
              read(99,*) i, c1(II,JJ)
              read(99,*) i, omega(II,JJ)
              read(99,*) i, b20(II,JJ)
              read(99,*) i,alpha2(II,JJ)

c     not used
c     read neutral paramters H3pO4
              read(99,*) i,xlamc
              read(99,*) i,xlama1
              read(99,*) i,xlama2
              read(99,*) i,xlamnn


              II=3
              jj=11

          filename=path(1:Lpath)//'/IVEA/model/b_Na3PO4.dat'           
           open(99,file=filename)
        II= 3
        JJ= 11
              read(99,*) i, b0(II,JJ)
              read(99,*) i, b1(II,JJ)
              read(99,*) i, c0(II,JJ)
              read(99,*) i, c1(II,JJ)
              read(99,*) i, omega(II,JJ)

c     set the 1 charged cations to Na
              DO II=1,nc
                 if (zc(ii).le.0.1) then
                  b0(II,JJ)= b0(3,JJ)
                  b1(II,JJ)= b1(3,JJ)
                  c0(II,JJ)= c0(3,JJ)
                  c0(II,JJ)= c0(3,JJ)
                  omega(II,JJ)= omega(3,JJ)
              endif
              enddo
              
        NC=6
        NA=11

        DO I=1,nc
        DO j=1,nA
              
              b0Z(i,j)=b0(i,j)
              b1Z(i,j)=b1(i,j)
              c0Z(i,j)=c0(i,j)
              c1Z(i,j)=c1(i,j)
              omegaZ(i,j)=omega(i,j)
           enddo
           
           enddo
           endif


        NC=6
        NA=11
        DO I=1,Nc
        DO j=1,NA
           omega(I,j)=omegaz(i,j)
           b0(I,j)=b0z(i,j)
           b1(I,j)=b1z(i,j)
           c0(I,j)=c0z(i,j)
           c1(I,j)=c1z(i,j)
        iflag(i,j)=0 
          
        enddo
        enddo
        


c     MA
c     1:NO3-
c     2:Cl-
c     3:HSO4-
c     4:SO4-2
c     5:CO3-2
c     6:HCO3-3

C     NC
c     1:H+
c     2:NH4+
c     3:Na+
        
        
c        nc=6
c        na=10
        
 
      	Iflag(1,1)=1            !H+ NO3-
	Iflag(1,2)=2 !H CL
	Iflag(2,1)=7            !NH4 NO3
	Iflag(2,2)=8  ! NH4 Cl
	NA=4
	NC=2
	Iflag(2,3)=5 ! NH4-HSO4
	Iflag(2,4)=6 ! NH4-SO4
	Iflag(1,3)=3 !H-HSO4
	Iflag(1,4)=4  ! H- SO4
      	call calpar(T,NC,NA,b0,b1,C0,C1,omega,xs,Iflag)

c     set  oxalatecto sulfate
c     set KHSO4 to NaHSO4


c        b0(3,5)= bb(1)
c        b1(3,5)= bb(2)
c        c0(3,5)= bb(3)
c        c1(3,5)= bb(4)
c        omega(3,5)= bb(5)
c        b20(3,5)=   bb(6)
c        alpha2(3,5)= bb(7)

        
        
	NA=11
	NC=6
        zc(4)=2d0
        zc(6)=2d0
c     replace NH4Cl
        b0(2,2)=0.29535E-02
        b1(2,2)=0.31808E+00  
        c0(2,2)=-0.14979E-04  
        c1(2,2)=0.49032E-01  
        omega(2,2)= 0.12116E+01 

        za(4)=2d0
        za(5)=2d0
        za(8)=2d0
        za(10)=2d0
        za(11)=3d0

        
c        DO I=1,Nc
c        DO j=1,NA
c           write(70,'(2I5,7E15.6)') i,j,b0(i,j), b1(i,j),c0(i,j),c1(i,j)
c     & ,omega(i,j),b20(i,j)          ,alpha2(i,J) 
c        enddo
c        enddo



        
        
        xm24=ml(24)
        xm29=ml(29)
        
        MC(1)= ml(6)
        MC(2)= ml(12)

c     treat  ML(19) as Na+ 
        MC(3)= ml(16) ! +ML(19)*1.0
         mc(4)= ml(32)
         MC(5)= ml(19)          !k
         MC(6)= ml(36)          !Ca++
         DO I=1,nc
            mc0(i)=mc(i)
            enddo
c     set the Na molality to maximal 40

        ma(5) =ml(14)
        ma(6) =ml(15)

        
        
        Ma(1)= ml(18)+ml(34)     !  treat lactic acid as HNO3
c     
        Ma(1)= ma(1) +ML(8)!  treat AA- as nitric acid

c     Treat ML(20) as Cl- 
        Ma(2)= ml(17)+ML(20)*1d0 

        ma(3)=ML(27)
        ma(4)=ML(28) !sum SO4--, HPO4--
        

        ma(7) =ml(25)
        ma(8) =ml(26)

        ma(9) =ml(30)
        ma(10) =ml(31)
        ma(11) =ml(37) + ml(38)

        

        do I=1,NA
        ma0(i)=ma(i)
        enddo
        xmm= 0
        DO I=1,NC
           xmm=xmm+mc(I)
        enddo
        
        DO I=1,NA
           xmm=xmm+ma(I)
        enddo
        
       DO I=1,NC
         Mmax=40d0
           if(I.eq.2) mmax=100

         if (mc(I).ge.mmax)  mc(I)=mmax
        enddo
        dd=3
        if (mc(4).ge.dd) mc(4)=dd
        if (mc(6).ge.dd) mc(6)=dd
        
        DO I=1,Na

         Mmax=40d0

         if(I.eq.1) mmax=100
           
           if (ma(I).ge.mmax)  ma(I)=mmax
           if (ma(8).ge.10) ma(8)=10d0
c           if (ma(26).ge.5) ma(26)=5d0
        enddo


        xmc= ml(15)+ml(16)
        ffc=1d0
        if (xmc.ge.15) ffc= 15/xmc
        ma(5) =ma(5)*ffc
        ma(6) =ma(6)*ffc

        

        
        if (xmm.le.1D-13) then
           awin=1d0
           gammah=1d0
           gammas1=1d0
           gammas2=1d0
           gammaNH4=1d0
           gammaNa=1d0
           gammacl=1d0
           gammaNO3=1d0
           gammaH3pO4=1d0
           gammaHOA=1d0
           gammaOA=1d0
           gammaCO3=1d0
           gammaHCO3=1d0
           gammaK=1
           gammaca=1
           gammamg=1
           gammaPO4=1
           
           goto 44
        endif
        
         
c        if (Ma(2).ge.Mmax) Ma(2)=Mmax
c        if (Ma(1).ge.Mmax) Ma(1)=Mmax
c        if (Ma(4).ge.35) Ma(4)=35

c     set the Cl molality to maximal 40
c        print*,'ma 4', ma(4), ma0(4)

        phinac= 0d0
        phiaa= 0d0
        
        
        etaNaCl=-0.002
        
      	awin=gammasn(T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega) ! Pitzer model ions
        xus=phiNaC*ma(5)*ma(6)*mc(3)+ ma(5)*ma(6)*phiaa !
c     CO2 mixture parameter
        
        xus=xus+ xlamna* mc(3)*mL(13)+
     & xlamcl* mA(2)*mL(13)+etanacl*ma(2)*mc(3)*ml(13)
        xme=0
        DO I=1,NC
           xme=xme+mc(I)
           enddo
        DO I=1,Na
           xme=xme+ma(I)
           enddo
        
        awin=awin*exp(-2/ml(1)*xus)

c        xmm=ma(1)+ma(2)+mc(1)+mc(2)
c        xmm0=ma0(1)+ma0(2)+mc0(1)+mc0(2)

        xmm=0
        xmm0=0
        DO I=1,2
        xmm=xmm+ma(I)
        xmm0=xmm0+ma0(I)
        enddo

        DO I=1,Nc
        xmm=xmm+mc(I)
        xmm0=xmm0+mc0(I)
        enddo


        awin=awin* (1000d0/MM(1)+ xmm)/(1000d0/MM(1)+xmm0)
 44     continue
        
        
c     calculates the aw of acetic acid, NH4CH4COO, and CO2

        xvH2O = 1000/MM(1)*MV(1)

c     vmol fraction of acetic acid and Co2 in ammonium acetate
c     and Phosphoric acid
        xmH2O=1000/mm(1)

c     all inorganic species not treated in Pitzer model
c        xvol = xmH2O/(XmH2O
c c    &       + ML(8)+ ML(9)+ ML(5) + ML(10) +ml (23) ) ! acetic acid ammonium acetate

        xvol=1d0

c     all neutral species
c        if (imode_output.eq.3) then
           
        xms=0


         xms=ml(2)
c       gammaco2= get_gammaco2(t,Ml)
        gammaH3Po4=dexp(xlama1*2 * ma(7) +2*xlama2*ma(8)
     &       + 2*xlamc *mc(1) + 2* xlamnn*ml(24))

         DO jj = 9, np-NPsolid
            dd=1
            if (jj.eq.13) dd= gammaco3
            if (jj.eq.24) dd= gammaH3po4
            
         if (Izc(jj).eq.0.and.(JJ.le.20 .or. jj.ge.23)) xms=xms+ ml(jj)
           enddo
c     special treatment of H3PO4 and organic acid 29
           xms=xms-ml(24)-ml(29)
           is_sucrose=0
           if (mm(2).ge.342.d0 .and. mm(2).le.343.4d0)is_sucrose=1
           if (is_sucrose.eq.1) then
              xms=xms-ml(2)
              endif
              aws = xmh2o/(xms+xmh2o) ! organics Raoult’s law

c     calculate the contribution of sucrose
              if (is_sucrose.eq.1) then
                 omega1= 1000/(1000d0+ mm(2)*ml(2))
                 call calaw_beni(T,omega1,awsuc)
                 aws =aws * awsuc
           endif
           




c     do empirical correction only for EDB
        if (imode_EDB.eq.1) then
        awin= aw_corr(awin)
        endif

        aw=awin*aws


        
c     print*, awin,aws,xvol
c        if  (imode_output.eq.0) then !ony valid for labo experiments
c        print*, 'b0', b0(NC,5), b0(NC,6)
c        print*, 'b1', b1(NC,5), b1(NC,6)
c        print*, 'c0', c0(NC,5), c0(NC,6)
c        print*, 'c1', c1(NC,5), c1(NC,6)
c        print*, 'o ', omega(NC,5), omega(NC,6)
c        print*, 'b20', b20(NC,5), b20(NC,6)
c        print*, 'a2 ', alpha2(NC,5), alpha2(NC,6)

c        DO I=1,NA
c           print*,i,mc(I),ma(i)
           
c        enddo
c        mc=0
c        ma=ma*20
c        mc=mc*20
c        DO xx=6.19,6.2,.3
c           ma=0
c           ma(5)= mc(nc)/2
        N1=1
	gammah=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)



        
        gammaH=gammaH*dexp(xlamc*2 * ml(24))
         
        N1=2
	gammaNH4=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        N1=3
	gammaNa=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
c        write(46,*) gammaNa
        xmix=phiNaC*ma(6)*ma(5)+etanacl*ml(17)*ml(13) 
        xmix=xmix+ 2 * xlamNA * ml(13)
                gammaNA=gammaNA*dexp(xmix)
        gamma2(Is,16) =gammaNa

        
        N1=5
	gammaK=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        gamma2(Is,19) =gammaK

        N1=4
	gammamg=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        gamma2(Is,32) =gammamg

        N1=6
	gammaca=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        gamma2(Is,36) =gammaca
        
        
        
        n2=NC+2
	gammacl=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        xmix=etanacl*ml(16)*ml(13)+2*xlamcl*ml(13)
        gammacl=gammacl*dexp(xmix)

        
        n2=NC+1
        gammaNO3=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

c     for pure phsospatem g=1

        n2=NC+3
        gammas1=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
c        ss= (ml(28)+ml(27)+ML(29)+ml(30)+ml(31))/(MA(3)+ma(4))

        
c     for pure phsospatem g=1
        n2=NC+4
        gammas2=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)


        

        n2=NC+5
        za(5)=2
        gammaCO3=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        n2=NC+6
        gammaHCO3=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        n2=NC+7
        gammaH2Po4=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        gammaH2Po4=gammaH2Po4*dexp(xlama1*2 * ml(24))

c        print*,xlamc,xlama1,xlama2,xlamnn

c        print*,'gammah3p', gammah3po4
        
        
        n2=NC+8
        gammaHPo4=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        gammaHPo4=gammaHPo4*dexp(xlama2*2 * ml(24))
        
        n2=NC+9
        gammaHOA=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        
        n2=NC+10
        gammaOA=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        n2=NC+11
        gammaPO4=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        
 
        ap = gammaNA*gammaCl*ml(16)*ml(17)
        snacl= ap/apnacl(T)
        xmm=0d0
        DO J=6,NP-npsolid
           xmm=xmm+ml(J)
        enddo
        xmm=xmm/2
         aa3= dsqrt(xmm)
         xlam_29=xlam1_29*dexp(-aa3/xmE1_29)+xlam2_29*dexp(-aa3/xme2_29)
         gammaH2OA = dexp(2* xlam_29*mL(29))

        xx= aa3/xmE1_29
	ggs=gs(Xx)
        xx= aa3/xmE2_29
	ggs2=gs(Xx)

        BSs =ggs*xlam1_29/xmm+ggs2*xlam2_29/xmm

            F2=Ml(29)*ml(29)*BSs
c            write(6,'(A,16E15.6)')'g29',ml(29),ml(2),gammaH2OA, dexp(F2)
c     &           , ggs,ggs2 ,xlam1_29, xlam2_29

            gammaH2OA =gammaH2OA *dexp(f2)
            
c     print*, ml(29), xlam_29, gammaH2OA

            if (ml(29).le.1D-30) gammaH2OA=1d0
            return
        end
      

c     ---------------------------------------------
c

c
c     HCl (gas )  --> H+(aq) + Cl-(aq)
c      
c     ---------------------------------------------

      function xhenry(T)
	
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 xhenry
	real*8 K0,K
	k0=2.05E6

	R=8.314
        T0=298.15
	dh=-74852.
	dcp=-165.52
	k=log(k0)+(dh-dcp*t0)/R*(1/T0-1/T)+ dcp/R*log(T/T0)
	k=exp(K)
	xhenry=k
	return
	end

c     ---------------------------------------------
c
c     The Henrys law coefficent of HNO3 in 1/atm, required 
c
c     HNO3 (gas )  --> H+(aq) + NO3-(aq)
c      
c     ---------------------------------------------
      
      
	function xkhx(T)

      IMPLICIT REAL*8(A-H,O-Z)
	real*8 xkhx
	x=385.9722-3020.3522/T-71.002*log(T)+.131442311*T
     * -.420928363E-4*T**2
	xkhx=exp(x)
	return
	end

c     -----------------------------------------------
c
c     calculate the activity of ions N (1 <= N <= Nc+Na) usin
c     Pitzer ion interactions model. 
c     
c     -----------------------------------------------
        function gammann(t,N,NC,NA,mC,mA,zC,zA,b0,b1,C0,C1,omega)
C	N: the index of the ion for which the activity is to be calculated
C	IF ( for cations: N <= Nc; else for anions N-NC )
C	zC, ZA: the ionic charge of type C or A.
C	b0,b1,c: Pitzer  coefficients B0(Nc, Na)
C	Nc : type of cations
C	Na : type of anions

C
C
C	This only a simple version, which only take the unsymmetrical 
C	factor e-theta and E-theta ' into account.
C	
C

	implicit real*8 (a-h,o-z)
	real*8 gammann
	Parameter( nmax=12 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 b(nmax,nmax),Bs(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
 	real*8 c(nmax,nmax),Cs(nmax,nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
	real*8 I, I2
	
	real*8 B20(nmax,nmax),alpha2(nmax,nmax)
	common /alpha/ alpha2, b20

	
c        print*,omega(1,1)

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
           if (mc(ic).gt.0d0) then
              I=I+mC(IC)*zC(IC)**2
              endif
	enddo
	DO IA=1,NA
           if (ma(ia).gt.0d0) then
	I=I+mA(IA)*zA(IA)**2
        endif
      enddo
	I=I/2.
        if (I.le.0d0) then
           gammann=1d0
           return
        endif

	I2=sqrt(I)

	z=0	
        gam=0.


	DO IC=1,NC
           if (mc(ic).gt.0d0) then
           
           Z=Z+ZC(IC)*MC(IC)
           endif
	enddo
	DO IA=1,NA
           if (ma(ia).gt.0d0) then

           Z=Z+ZA(IA)*MA(IA)
           endif
	enddo



CCC	Calculate B,Bs, C, Cs
	x=sqrt(I)*alpha
	gg=g(x)
	ggs=gs(X)

	DO Ic = 1, NC
	DO Ia = 1, Na
           if (mc(ic).gt.0d0 .and. ma(ia).gt.0d0) then
           if (alpha2(Ic,ia).le.0d0) alpha2(Ic,ia)=1d0
		x2=sqrt(I)*alpha2(ic,ia)
		gg2=g(x2)
	        ggs2=gs(X2)
	B(IC,IA) = b0(Ic,Ia) + gg* b1(Ic,Ia)+ gg2*b20(Ic,iA)
   	BS(IC,IA) = ggs*b1(Ic,Ia)/I+ ggs2*b20(Ic,Ia)/I

           
c         B(IC,IA) = b0(Ic,Ia) + gg* b1(Ic,Ia)
c	BS(IC,IA) = ggs*b1(Ic,Ia)/I
	omega1=omega(Ic,IA)
        xo=omega1*sqrt(I)
	xhx=1./xo**4*(6.-exp(-xo)*(6+6*xo+3*xo**2+xo**3))
	xhxs=exp(-xo)/2-2*xhx

	C(Ic,Ia)=(C0(Ic,Ia)+4*C1(Ic,Ia)*xhx)
	Cs(Ic,Ia)= C1(Ic,Ia)/I*xhxs
        endif
      enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2

	F1= -Aphi*(I2/(1.+1.2*I2) +2./1.2*log(1+1.2*I2))
	F2=0

	DO Ic = 1, NC
	DO Ia = 1, Na
           if (mc(ic).gt.0d0 .and. ma(ia).gt.0d0) then


              F2=F2+mc(Ic)*Ma(Ia)*(BS(IC,Ia)+2*Z*CS(Ic,IA) )
              endif
           enddo		
	enddo

	xm=mc(1)
c	J1=1
c	J2=2
c	call EFUNC(J1,J2,Aphi,I,E,ED)
	
	F3= 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
c        z1=ZC(IC1)c
c     z2=ZC(IC2)
           if (mc(ic1).gt.0d0 .and. mc(ic2).gt.0d0) then
           
           J1=ZC(IC1)+.1
           J2=ZC(IC2)+.1
	IF(j1.eq.j2) goto 21
	call EFUNC(J1,J2,Aphi,I,E,ED)
	F3=F3+ED*MC(Ic1)*MC(IC2)
21	continue
        endif
      enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
c	z1=ZA(IA1)
c	z2=ZA(IA2)
           if (ma(ia1).gt.0d0 .and. ma(ia2).gt.0d0) then

           J1=Za(Ia1)+.1
           J2=Za(Ia2)+.1
	IF(j1.eq.j2) goto 22
	call EFUNC(J1,J2,Aphi,I,E,ED)
	F4=F4+ED*MA(IA1)*MA(IA2)
22	continue
        endif
      enddo
	enddo
	

	F=F1+F2+F3+f4

c	print*, 'F=', F
	


	IF (N.le.NC) then

	F3=0
	DO Ic1 = 1, NC
c	z1=ZC(N)
c     z2=ZC(IC1)
           if (mc(ic1).gt.0d0 ) then
           
           J1=ZC(N)+.1
           J2=ZC(IC1)+.1
	IF(j1.eq.j2) goto 121
	call EFUNC(J1,J2,Aphi,I,E,ED)

	F3=F3+E*MC(Ic1)
 121    continue
        endif
	enddo

	a1=zc(N)**2 *F
	a2=0
	DO Ia=1,NA
           if (ma(ia).gt.0d0) then
           a2=a2+MA(IA)*( 2*B(N,IA)+Z*C(N,IA) )
           endif
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
           if (mc(ic).gt.0d0 .and. ma(ia).gt.0d0) then
                      a3=a3+MA(IA)*MC(IC)*C(IC,IA)
                      endif
                   enddo
	enddo

	a3=Zc(N)*A3
	gam=a1+a2+a3 +F3


	else 

	N1=N-Nc
	a1=zA(N1)**2 *F

	f4=0.
	DO IA1 = 1, NA
c	z1=ZA(N1)
c     z2=ZA(IA1)
           if ( ma(ia1).gt.0d0) then
           
           J1=ZA(N1)+.1
           J2=Za(IA1)+.1
	IF(j1.eq.j2) goto 122
	call EFUNC(J1,J2,Aphi,I,E,ED)
	F4=F4+E*MA(IA1)
122	continue
        endif
      enddo

	
	a2=0

	DO IC=1,NC
           if (mc(ic).gt.0d0 ) then

           a2=a2+MC(IC)*( 2*B(IC,N1)+Z*C(IC,N1) )
           endif
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
           if (mc(ic).gt.0d0 .and. ma(ia).gt.0d0) then
           a3=a3+MA(IA)*MC(IC)*C(IC,IA)
           endif
	enddo
	enddo

	a3=ZA(N1)*A3
	gam=a1+a2+a3+f4
c	print*, 'gam', n1,a1,a2,a3,f4,za(n1)

      endif
	gammann= exp(gam)
	return
	end
        

c     -----------------------------------------------
c
c     calculate water activity of ionic solution using
c     Pitzer ion interactions model. 
c     
c     -----------------------------------------------
      


      
c     ----------------------------------------------------------
c      required by the Pitzer model gammann and gammasn
c     ----------------------------------------------------------
      function gs(X)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 gs
	gs=2*(-1.+(1+x+x**2/2.)*exp(-x) )/x**2
	return
	end


c     ----------------------------------------------------------
c      Liniear interpolation subroutine
c     ----------------------------------------------------------

      subroutine intpl(x1,y1,n1,x2,y2,n2)
           implicit real*8 (a-h,o-z)
           REAL*8 X1(*),X2(*),y1(*),y2(*)
           y2(1)=y1(1)
           
           DO I=1,n2

              if(x2(i).le.x1(1)) then
                 y2(i)=y1(1)
                 goto 20
              endif

              if(x2(i).ge.x1(n1)) then
                 y2(i)=y1(n1)
                 goto 20
              endif


              DO J=2,n1
                  if(x1(J).ge.x2(I)) then
                     ff= (x2(I)-X1(J-1))/(x1(J)-x1(J-1))
                  yy=y1(J-1)+ff*(y1(J)-y1(J-1))
                  y2(I)=yy
                  goto 20
                  ENDIF

             enddo

 20           continue

           enddo



           return
           end

c     ----------------------------------------------------------
c      subroutine to find the zero place of function func
c     ----------------------------------------------------------
c     this code works fun(a) > 0
c     this code works fun(b) < 0




c     ----------------------------------------------------------
c      subroutine to find the zero place of function func, for the case when they used twice
c     ----------------------------------------------------------
c     ---------------------------------------------
c      
c     Generate the parameters of ion paars for Pitzer ion interaction model
c
c     
c     --------------------------------------------
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this program is valid only for the interaction, where at least C
C     one ion is single charged.                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



  
C     Iflag              Molality range          T-range
C     1: H-NO3            0-16                   190-320
C     2: H-Cl             0-60                   190-320
C     3: H-HSO4           0-40                   190-320 
C     4: H-SO4            0-40                   190-320
C     5: NH4-HSO4         0-30                   220-320
C     6: NH4-SO4          0-30                   220-320
C     7: NH4-NO3          0-10                   220-320
C     8: NH4-Cl           0-7                    220-320 
C     9: Na-SO4		  0-4                    220-320
  


C	xs is the mixture parameter
C	1: psi(H,HSO4,NO3)	
C	2: theta(SO4,NO3)	
C 	3: psi(H,SO4,NO3)
C	4: psi(H,HSO4,Cl)	
C	5: theta(HSO4,Cl)	
C 	6: psi(H,SO4,Cl)
CC	7: psi-NH4-SO4-HSO4
CC	8: psi-NH4-H-SO4
CC	9: psi-NH4-H-HSO4
CC     10: Phi-NH4-H
CC     11: psi NH4-NO3-SO4 
CC    


   	subroutine calpar(T,NC,NA,b0,b1,C0,C1,omega,xs,Iflag)
  
	Parameter( nmax=12 )
        IMPLICIT REAL*8(A-H,O-Z)
  	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
	integer Iflag(NMAX,NMAX)
 	real*8 bb(50),b2(50),b3(50),b8(50)
	real*8 b5(50),b6(50),b7(50),omega(nmax,nmax),xs(100),b9(50)

CCCCCCCCCCCC data for H-Cl CCCCCCCCCCCCCCCCCCCC
        data (b2(I),i=1,21) /
     > 0.23378,-7.21238E-02,-1.7335667E-02,5.760665E-03,-8.29279E-03,
     >0.2897,7.575434E-02,-1.1474E-03,0.38038,-0.309442,-2.794885E-03,
     >2.309349E-04,9.322982E-04,-2.398E-04,2.85959E-04,-0.21154,
     > 0.101481,5.945618E-02,-0.107864, 8.81749E-02,  1.9916/
       
CCCCCCCCCCCCC data for H-NO3 CCCCCCCCCCCCCCCCC
        data (bb(I),I=1,27)/
     * 3.895835E-03,-1.55571E-02,1.703729E-02,
     * -5.6173712E-03,  5.732047E-03,  0.91622,  0.613523,
     * -0.68489, 0.3038,  -0.32888, 7.6086113E-07, 7.2714678E-05,
     * -1.0037E-04,3.475E-05,-3.62927E-05,5.380465E-02,-2.2163E-02,
     * -1.0166E-02, 6.5423E-03,  -8.80248E-03, 0.907342,
     *         -6.78428E-4,9.576E-4,2*0D0,  7.769E-3, -5.819E-4 /

c      !mixturedata
   
CCCCCCCCCCCC data for H-HSO4, H-SO4 CCCCCCCCCCCCCCCCCCCC
        data (b3(I),i=1,42) /  0.148843, -7.769E-2 ,  2.8062E-2 ,
     *  4.7903E-4,  7.25E-4 ,  0.17843,  0.678,  8.7381E-2, 
     *  -0.57881,  7.58E-2 , -9.878E-4 ,  5.447651E-4, -2.58798E-4,
     *   1.8466527E-5 ,  1.23457E-5 ,  0.37138, -9.24874E-2 ,
     *  -9.21372E-3 , -1.065158E-2 ,  5.4987733E-2 ,  0.2726312,
     * -1.34824E-3 , -0.24711,  1.25978E-2 ,  0.11919,  0.7397,
     *  -3.01755,  -4.5305,  -3.1072, -0.8555842,  9.2223E-4 ,
     * -4.1694532E-3 ,  7.141266E-3 ,  2.32984E-3 , -6.98191E-4 ,
     *  -2.242,  0.71925,   2.52, -0.7391,  -1.548503, 1.5452, 2./


CCCCCCCCCCCCCCCC NH4-HSO4 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     with out chan's data

c	data (b5(I),I=1,32)/
c     * -8.746E-4,  -2.3125, -9.56785E-6,   2.58238,   2.38, -3.1314E-4,
c     *  1.6896E-2, -0.7351,  0.6883,  1.813E-3, -0.1012515 ,
c     * -2.66E-2, -2.86617E-3,  0.22925,  0.438188,  2.522E-4,
c     * -2.90117E-5,  0.9014,  0.41774,  -1035.9,  0.0,  -299.69,
c     *  0.0, -4.9687E-4,  0.0,  1.21485E-2,  0.0, -1.0334E-3,
c     *  0.0,  8.48374E-2,  0.0,0d0/

	data (b5(I),I=1,32)/
     &    -0.11671E-01,  
     2   -0.49047E+00 ,  
     3    0.10378E-03 ,  
     4    0.10170E+00 ,  
     5    0.11722E+01 ,  
     6   -0.30795E-03 ,  
     7   -0.20000E-03 ,  
     8   -0.35877E+01 ,  
     9    0.23663E+01 ,  
     1   -0.10000E-03 ,  
     1    0.50000E-04 ,
     &	  0.169702108786968   
     *, -0.387640333038779    
     *,  0.0
     *,  0.0
     *,  2.601561402128945E-004
     *, -2.337658723484350E-003
     *,  0.0
     *,  0.0
     *,  -612.096813408415 
     *,  0.00
     *,  -306.048406704207     
     *,  0.00
     *,  5.320349188422027E-004
     *,  0.0
     *,  0.0
     *,  0.0
     *,  3.532241943713784E-003 
     *,  0.0
     *, -4.534300563622422E-002 
     *,  0.0000
     *,  0.0000/
	data ikk /0/
	save ikk
	
	ikk=1

      
C      ! with chan's data
c	data (b5(I),I=1,31)/
c     * -7.8224E-3,  -1.722,  7.5882E-5,  0.962,   1.8914, -4.285E-4,
c     *  3.99E-4, -0.76151,  0.53133,  1.16726E-3,
c     * -9.45855E-2, -2.007E-2,  6.0482E-3, -0.1392,
c     *   2.9544,  2.45133E-4,   -5.492E-5,  0.8358,
c     * -0.3791,  -980.81,  0.,   265.7,  0., -1.7568E-3,
c     * -7.72E-4,  3.437E-4,   -1.6456E-2, -1.2034E-3,
c     * -2.2943E-3,  7.9281E-2,  3.772E-2/

CCCCCCCCCCCCCCCC NH4-SO4 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	data (b6(I),I=1,13) /-1.2058223E-2 ,1.1043,4.79018E-5
c     * ,2.14346E-2 , 0.58,  -2.9146E-2,  1.9631E-4  ,   1.1378,
c     *  0.9283,1.28548E-4,1.684E-5,2.6267E-2, -2.6E-4/   ! wt=1
 
C	data (b6(I),I=1,13)/
C     * -1.2058223E-2,   1.1043,  4.79018E-5,2.14346E-2 ,  0.58,
C     * -0.1188,  8.5E-2,   2.10514,  0.5942,  7.888E-4, -5.503E-4,
C     *  5.815E-2 , -3.766E-2 / !wt=5d-5

	data (b6(I),I=1,13)/
     *  -4.327681689677379E-003 
     *      ,  0.953787255317688      
     *      ,  1.113167041454729E-005 
     *      ,  2.101340928698026E-002 
     *      ,  0.614050552546192      
     *      , -2.606487385557435E-003 
     *      ,  5.672030237892362E-004 
     *      ,   1.75633200902706      
     *      ,   1.10907828739749      
     *      ,  1.177526981905975E-005 
     *      ,  6.713875789826932E-007 
     *      ,  5.973238697090055E-003 
     *      , -2.812698651271791E-003/ 
	 



CCCCCCCCCCCCC data for NH4-NO3 CCCCCCCCCCCC
	data (b7(I),I=1,13) /-2.3275E-2, 0.15, 1.1634E-4,1.62E-3,
     * 0.43, ! 0.107264,  0.221, -3.8442E-4, -1.3872E-2, 4*0D0/
     *  8.78E-2 ,  0.2753645, -3.349E-4 , -1.093E-2 ,
     * -4.769E-2 ,  0.1776,  1.25E-4 ,  6.9751E-3 /

C
C     *  2.15438E-2,0.67073,  2*0.0, -1.3662E-2,3.4747E-002,2*0.0/
 

CCCCCCCCCCCCCCCC NH4-Cl CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	data (b8(I),I=1,9) / -6.333E-4 , -3.99546E-4 , 0.3155,
     *  0.1414, -3.837E-5, 1.08331E-4, 5.2436E-2,1.6827E-2,1.19/

CCCCCCCCCCCCC data for Na-SO4 CCCCCCCCCCCC
	data (b9(I),I=1,17) /1.63E-03 , 2.092E-03,3.484156E-02 ,
     * -1.057E-02,1.0775,  0.9, 0.8206,-9.6425E-02,2.7492E-03,
     *-9.2838E-04,-6.8268E-04,4.8126E-05,7.182E-02,
     *  0.7586, -0.30291, 0.2311,1.7/

	xs(4)=0
	xs(5)=0
	xs(6)=0
	DT=(T-298.15)/100D0

	DO I=1,NC
	DO J=1,NA
	if (Iflag(I,J) .eq. 2) then

	omega(i,j)=b2(21)
        B0(i,j)=b2(1)+dt*b2(2)+dt**2*b2(3)+dt**3*b2(4)+dt**4*b2(5)
        B1(i,j)=b2(6)+dt*b2(7)+dt**2*b2(8)+dt**3*b2(9)+dt**4*b2(10)
        c0(i,j)=b2(11)+dt*b2(12)+dt**2*b2(13)+dt**3*b2(14)+dt**4*b2(15)
        c1(i,j)=b2(16)+dt*b2(17)+dt**2*b2(18)+dt**3*b2(19)+dt**4*b2(20)
        endif

	 IF (Iflag(I,j).eq.7) then
	b0(I,J)=b7(1)+b7(6)*dt+b7(10)*dt*dt
	b1(I,J)=b7(2)+b7(7)*dt+b7(11)*dt*dt
	c0(I,J)=b7(3)+b7(8)*dt+b7(12)*dt*dt
	c1(I,J)=b7(4)+b7(9)*dt+b7(13)*dt*dt
	omega(I,j)=b7(5)

	xs(1)=bb(22)+bb(23)*dt
	xs(2)=bb(24)+bb(25)*dt
	xs(3)=bb(26)+bb(27)*dt
        xs(11)=4.75458E-4-4.0577E-003*dt
c4.633E-4-4.093E-3*dt
        endif

	IF (Iflag(I,J).eq.9) then

        B0(i,j)=b9(1)+dt*b9(2)+dt**2*b9(3)+dt**3*b9(4)
        B1(i,j)=b9( 5)+dt*b9(6)+dt**2*b9(7)+dt**3*b9(8)
        c0(i,j)=b9( 9)+dt*b9(10)+dt**2*b9(11)+dt**3*b9(12)
        c1(i,j)=b9(13)+dt*b9(14)+dt**2*b9(15)+dt**3*b9(16)
	omega(I,j)=b9(17)
        endif

	IF (Iflag(I,j).eq.8) then

	omega(I,j)=b8(9)
	b0(I,J)=b8(1)+b8(2)*dt
	b1(I,J)=b8(3)+b8(4)*dt
	c0(I,J)=b8(5)+b8(6)*dt
	c1(I,J)=b8(7)+b8(8)*dt
        endif

	 IF (Iflag(I,j).eq.3) then
	omega(I,J+1)=b3(42)
	dt=(t-298.15)/100.
	omega(I,J)=b3(41)
        B0(i,j)=b3(1)+dt*b3(2)+dt**2*b3(3)+dt**3*b3(4)+dt**4*b3(5)
        B1(i,j)=b3(6)+dt*b3(7)+dt**2*b3(8)+dt**3*b3(9)+dt**4*b3(10)
        c0(i,j)=b3(11)+dt*b3(12)+dt**2*b3(13)+dt**3*b3(14)+dt**4*b3(15)
        c1(i,j)=b3(16)+dt*b3(17)+dt**2*b3(18)+dt**3*b3(19)+dt**4*b3(20)
      endif
      
	 IF (Iflag(I,j).eq.4) then
	omega(I,J)=b3(42)
       B0(i,j)=b3(21)+dt*b3(22)+dt**2*b3(23)+dt**3*b3(24)+dt**4*b3(25)
       B1(i,j)=b3(26)+dt*b3(27)+dt**2*b3(28)+dt**3*b3(29)+dt**4*b3(30)
       c0(i,j)=b3(31)+dt*b3(32)+dt**2*b3(33)+dt**3*b3(34)+dt**4*b3(35)
       c1(i,j)=b3(36)+dt*b3(37)+dt**2*b3(38)+dt**3*b3(39)+dt**4*b3(40)
        endif

	if ( Iflag(I,J) .eq. 1) then
	omega(i,j)=bb(21)
	dt=(T-298.15)/100.
        B0(i,j)=bb(1)+dt*bb(2)+dt**2*bb(3)+dt**3*bb(4)+dt**4*bb(5)
        B1(i,j)=bb(6)+dt*bb(7)+dt**2*bb(8)+dt**3*bb(9)+dt**4*bb(10)
        c0(i,j)=bb(11)+dt*bb(12)+dt**2*bb(13)+dt**3*bb(14)+dt**4*bb(15)
        c1(i,j)=bb(16)+dt*bb(17)+dt**2*bb(18)+dt**3*bb(19)+dt**4*bb(20)


	xs(1)=bb(22)+bb(23)*dt
	xs(2)=bb(24)+bb(25)*dt
	xs(3)=bb(26)+bb(27)*dt
          xs(11)=4.75458E-4-4.0577E-003*dt
c4.633E-4-4.093E-3*dt
        endif

	if ( Iflag(I,J) .eq. 5) then

	B0(i,j)=b5(1)+b5(11+1)*(dt)+b5(11+2)*dlog(t/298.15)
	B1(i,j)=b5(2)+b5(11+3)*(dt)+b5(11+4)*dlog(t/298.15)
	c0(i,j)=b5(3)+b5(11+5)*(dt)+b5(11+6)*dlog(t/298.15)
	c1(i,j)=b5(4)+b5(11+7)*(dt)+b5(11+8)*dlog(t/298.15)
	omega(i,j)=b5(5)+b5(32)*dt

	xs(7)=b5(6)  +b5(11+13)*(dt)+b5(11+14) *dlog(t/298.15)
	xs(8)=b5(7)  +b5(11+15)*(dt)+b5(11+16) *dlog(t/298.15)
	xs(9)=b5(10) +b5(11+17)*(dt)+b5(11+18) *dlog(t/298.15)
	xs(10)=b5(11)+b5(11+19)*(dt)+b5(11+20) *dlog(t/298.15)
c	xs(10)=0d0

	T0=298.15

        sld1 = -86.+2791.9/T+13.482*dlog(T)
	sld2=b5(8)+(1/t-1/298.15)*b5(11+9)+b5(11+10)*log(t/t0)
	sld3=b5(9)+(1/t-1/298.15)*b5(11+11)+b5(11+12)*log(t/t0)
        endif

	if ( Iflag(I,J) .eq. 6) then

	omega(I,J)=b6(5)
	b0(I,J)=b6(1)+b6(6)*dt+b6(7)*dt*dt
	b1(I,J)=b6(2)+b6(8)*dt+b6(9)*dt*dt
	C0(I,J)=b6(3)+b6(10)*dt+b6(11)*dt*dt
	C1(I,J)=b6(4)+b6(12)*dt+b6(13)*dt*dt
                   dd=00

	endif


	enddo
	enddo

	xs(11)=0d0
	
	return
	end

c     ---------------------------------
c
c     vapor pressure of pure water Koop and Murphy
c
c     -----------------------------------
      
      function vwater(temp)
      implicit real*8 (a-h,o-z)
      vwater= dexp(54.842763 - 6763.22/temp - 4.210*dlog(temp) +
     + 0.000367*temp + dtanh(0.0415*(temp - 218.8))*(53.878 -
     + 1331.22/temp - 9.44523*dlog(temp) + 0.014025*temp))
	vwater=vwater/100d0

	return
        End

c     ---------------------------------
c
c     required by Pitzer model gammann and gammasn
c
c     -----------------------------------

	function g(x)
      IMPLICIT REAL*8(A-H,O-Z)
      		real*8 g
      	g=2*(1.-(1.+x)*exp(-x))/X**2
	return
	end

c     ---------------------------------
c
c     required by Pitzer model gammann and gammasn
c
c     -----------------------------------

	SUBROUTINE EFUNC(Icharg,JCHARG,A,XI,E,ED)
       IMPLICIT REAL*8(A-H,O-Z)
		REAL*8	X(3),J0(3),J1(3),DUM

	IF((ICHARG.EQ.JCHARG) .OR. (XI .LE. 1.D-30) )THEN
	E=0
	ED=0
		ELSE

	X(1)=6*ICHARG*JCHARG*SQRT(XI)*A
	X(2)=6*ICHARG*ICHARG*SQRT(XI)*A
	X(3)=6*JCHARG*JCHARG*SQRT(XI)*A

	DO I=1,3		
	DUM=-1.2D-2*X(I)**.528
	J0(I)=X(I)/(4.+4.581*X(I)**(-.7238)*EXP(DUM))
	J1(I)=(4.+4.581*X(I)**(-.7238)*EXP(DUM)*(1.7238-DUM*.528
     * ))/(4.+4.581*X(I)**(-.7238)*EXP(DUM))**2
	ENDDO
	
	XE=ICHARG*JCHARG/4./XI*(J0(1)-.5*J0(2)-.5*J0(3))
	XED=ICHARG*JCHARG/8./XI**2*(X(1)*J1(1)-.5*X(2)*J1(2)-.5*
     * X(3)*J1(3))-XE/XI
	E=XE
	ED=XED
	ENDIF
	
	RETURN
	END

c     --------------------------------------------------------
c
c      Reset shells to equal volume when one shell becomes too thin or too thick
c
c     -------------------------------------------------------
      subroutine reset_shells(NS,x)
       implicit real*8 (a-h,k,m, o-z)
      parameter (NSMM=100)
      parameter (np=50, npsolid=12)         ! number of species
      real*8 xn(NSMM,np),xn0(NSMM,np),vshell(NSMM)
      common /xn/xn0


      real*8 xn2(NSMM,np),ml(np)
      real*8 x(*), x0(NSmm+1)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      
	common /m/ mm, mv,izc(NP)
      common/solid/xnsolidt

c      common /vshell/ vshell
      common /pi/ pi
      
c      xn=xn0
      
      Do I=1,NS
         DO J=1,NP
      xn(i,j)=xn0(i,j)
      enddo
      enddo
      
c      x0(1)=0d0

c     set equidstance shells to equil volume

c      print*, 'diff tau = ', difftau
      volmin=1D8
      volmax=0
      
      DO I=1,NS
         vv = mv(1)*xn0(I,1)+mv(2)*xn0(I,2)

          DO j=6,np-NPsolid
          if (j.lt.21 .or. j.gt. 22)  
     & vv = vv +mv(j)*xn0(I,j)
             enddo
             
         if (vv.gt. volmax) volmax= vv
         if (vv.lt. volmin) volmin= vv
         
         enddo
         ff=1.5

         if (xnsolidt.ge.1D-24) ff=2
         if (volmax.gt. ff* volmin .and. volmax.gt.0d0) then
            print*, ' Rebin '
            write(28,'(A)') ' rebin'


       DO I=1,ns
          DO J=1,NP-npsolid
          xn(I,j)=0d0
          enddo
       enddo


       vcore=4*pi/3*x(1)**3
       vol=vcore
             x0(1) = (vol/4d0/pi*3d0 )**(1/3d0)
       DO I=1,ns
          vol = vol+mv(1)*xn0(I,1)+mv(2)*xn0(I,2)
          vv = mv(1)*xn0(I,1)+mv(2)*xn0(I,2)

          DO j=6,np-npsolid

          if (j.lt.21 .or. j.gt. 22)  
     & vol = vol +mv(j)*xn0(I,j)
          if (j.lt.21 .or. j.gt. 22)  
     & vv = vv +mv(j)*xn0(I,j)

             enddo
             x0(I+1) = (vol/4d0/pi*3d0 )**(1/3d0)
c             print*,i,x0(I+1), x(I+1)
c             print*,vshell(I), vv

             enddo

             vol1=(vol-vcore)/NS

             print*,' vol1 = ', vol1

     

       DO I=1,NS
C     CALCULATE XN1,XN2
          X(I+1)= ((vcore+vol1*I)/4d0/pi*3d0)**(1/3d0)
          
          RI = X(I)
          RA = X(I+1)

          DO J=1 ,ns
             if (x0(J+1).ge. RI) goto 11
          enddo
 11       NI= J

          DO J=1 ,ns
             if (x0(J+1).ge. RA) goto 12
          enddo
 12       NA= J
           if (NA.gt.NS)NA=NS

c          print*,I, NI,NA,x(I+1)



          if (NA.eq.NI)  then
          vnew = (RA**3  - RI**3)
          vS=( x0(NI+1)**3 -x0(NI)**3 )
          DO J=1,NP-npsolid
          xn(I,j)= vnew/vs* xn0(NI,j)
          enddo
          goto 13
          endif


c     take the RI to x(NI+1)
          vnew= x0(NI+1)**3 -rI**3 
          vS= x0(NI+1)**3 -x0(NI)**3 

          DO j=1,Np-npsolid
          xn(I,J)=xn0(NI,j)*vnew/vs
          enddo

c     take x(NA) to RA
          vnew=( RA**3 -x0(NA)**3 )
          vS=( x0(Na+1)**3 -x0(Na)**3 )

          DO j=1,Np-npsolid
          xn(I,j)=xn(I,j)+ xn0(NA,j)*vnew/vs
          enddo

C     TAKE THE SHELLS NI+1 TO NA-1npsolid
          DO jJ = ni+1,NA-1
             DO j=1,np-npsolid
          xn(I,j) =xn(I,j)+xn0(jj,j)
          enddo
          enddo
13       continue
       ENDDO



       DO J=1,NP
c         write(28,'(I5,100E15.6)') J, (xn(I,j),I=1,ns),(xn0(I,j),I=1,ns)

          xns0=0d0
          xns=0d0
         DO I=1,ns
            xns=xns+ xn(i,j)
            xns0=xns0+ xn0(i,j)
          enddo
c          if (xns0.gt.0d0) print*,j, xns, (xns-xns0)/xns0
c          if (xns0.gt.0d0) write(28,*) j, xns, (xns-xns0)/xns0

          enddo
          print*, 'finish rebin'          


         vol=0d0
         vol0=0d0

c     Reset  the bins for virus
          DO I=1,NS
             xn(i,21) = xn0(i,21) 
             xn(i,22) = xn0(i,22) 
          enddo

            vol=0
            DO I=1,NS
               DO J=1,np
            xn0(i,j)=xn(i,j)
         enddo
         enddo
         
            T=298.15
            ML(1)=1000D0/mm(1)
            DO I=1, NS
               DO J=1,np
                  ML(J)= xn(I,j)*ml(1)/xn(I,1)
                 
               enddo
                  
           call calHNew(T,ml)
                      call aw_back
     & (t,Ml,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
                      phh= dlog(mL(6)* gammah)/dlog(.1d0)
                      print*, 'phh = ', i,phh

               enddo
          DO I=1,NS+1
c             print*,i, x(I)
             enddo
c          stop

            endif
c
          
          
       return
       end
 

c     --------------------------------------------------------
c
c      Merge shell when it becomes very thin
c
c     -------------------------------------------------------

      subroutine merge(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np),xna(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),ml(np),flux(NSMM,np+npsolid)
      
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell
     & ,vshell0(nsmm)
c     common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa

      real*8 dx1(nsmm)
      
       common /dxmin/dxmin0
      common /NS0/NS0
       if (NS.le.1) return

      vmax=0
      DO I=1,NS
         vol(i)= xn(I,1)*mv(1)+xn(I,2)*mv(2)
         DO J=6,np-npsolid
         vol(I)=vol(I)+xn(I,j)*mv(j)
         enddo
      enddo

      vmax=0
      
      dxmin=dxmin0*.75
c      print*, dxmin, dxmin0
      NS1=ns
      DO JJ=1,ns1
c      if (NS.le.ns0) return

      rcore3= x(1)**3
      DO I=ns,1,-1
         
         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( rcore3+vs/4/pi*3)**(1/3d0)

         vdry= vs+ xn(I,2)*mv(2)+mv(1)*xn(I,1)

       DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rr= ( rcore3+vdry/4/pi*3)**(1/3d0)
         dx=rr-xss
         dx1(I)=dx
      enddo
      
c         print*,' dx, dxmin', I,dx, dxmin
      dminn=11110d0
      DO I=1,NS
         if (dx1(I).le.dminn) then
            dminn = dx1(I)
            Imin=I
         endif
         
      enddo
      
      DO IiI =1,1
      I = Imin


         

         if (dx1(I).le.dxmin) then

            I1= I


            if (I.gt.1.and.I.lt.NS) then
               if (vol(I-1).ge.vol(I+1)) I1=I-1
            endif
            if (I1.ge.ns) I1=Ns-1
            
            
            if(NS.le.ns0) then
           print*, 'merge  redistribute ', I1,dx
           xna0=0
           DO II=1,NS
              xna0=xna0+xn(II,16)
           enddo
           DO J=1,np-npsolid
              ss=xn(I1,j)+xn(I1+1,j)
              xn(I1,j)=ss/2
              xn(I1+1,j)=ss/2
           enddo
           write(6,'(A,100E15.6)')'merge redis',(dx1(kk),kk=1,ns)
            goto 222
               
            endif
            
            
           print*, 'merge dx ', ns,dx
           
            
            DO J=1,np
               xn(I1,J)=xn(I1,J)+xn(I1+1,J)
            enddo
            X(i1+1)= X(i1+2)

            ns=ns-1
            do iI=i1+1,NS
               do j=1,nP
               xn(II,J)=xn(II+1,J)
            ENDDO
            X(ii+1)= X(ii+2)
            ENDDO
            
            goto 222
         endif
         
c     5*dx
         dd=x(NS+1)-x(1)/10
c         print*,'ns, dd,dx,', ns,dd,dx
         
         if (NS.ge.10) then
            if (dx.le. dd) then
             if (I.gt.1) then

                dh=dabs( phshell(I)-phshell(I-1))
               if (i.lt.NS) dh1=dabs( phshell(I)-phshell(I+1))
               dh1=0
               if (dh1.gt. dh) dh= dh1
            else
              dh=dabs( phshell(I)-phshell(I+1))
               
           endif
           
             if (I.gt.1) then
                daw=dabs( awshell(I)-awshell(I-1))
                daw1=0d0
                if (i.lt.NS) daw1=dabs( awshell(I)-awshell(I+1))

               if (daw1.gt. daw) daw= daw1
            else
              daw=dabs( awshell(I)-awshell(I+1))
               
            endif

            
c
            print*,'dh, daw ', I,dh,daw
            
            if (dh.le. 0.05 .and.daw.le..01) then
               print*, 'merge, 5Xdxmin0m', I, dx, dh
            I1= I
            if (I.eq.NS) I1=I-1
            if (I.eq.1) I1=1
            if (I.gt.1.and.I.lt.NS) then
               I1= I
               if (vol(I-1).lt.vol(I+1)) I1=I-1
            endif
            

            DO J=1,np
               xn(I1,J)=xn(I1,J)+xn(I1+1,J)
            enddo
            X(i1+1)= X(i1+2)

            ns=ns-1
            do iI=i1+1,NS
               do j=1,nP
               xn(II,J)=xn(II+1,J)
            ENDDO
            X(ii+1)= X(ii+2)
            ENDDO
            
            goto 222


            endif
            

            endif
         endif
         

      enddo
      return
 222  continue
                     xnca=0
             DO II=1,NS
         xnca=xnca+ xn(Ii,16)
      enddo
      
         npl=np-npsolid
c                     timea=1
c                     call cal_flux(timea,x,flux,dtime,NS)
                     
             DO II=1,NS


                ml(1)=1000/mm(1)
                
               DO J=2,np
                  ml(j)= ml(1)*xn(II,J)/xn(II,1)

               enddo
c               DO KK=1,10
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

        
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
          awshell(II)= aw
c          print*,ii, aw,x(II+1)-x(II)
       enddo

c            enddo
            print*,'after merge d x na' ,xnca,xna0
            return

         enddo

      if (NS.le.ns0) return
         
         DO I=1,NS
c     compare shell below 
      ismerge =0

             if (I.gt.1) then
               dh=dabs( phshell(I)-phshell(I-1))
               if (i.lt.NS) dh1=dabs( phshell(I)-phshell(I+1))
               if (dh1.gt. dh) dh= dh1
               daw=dabs( awshell(I)-awshell(I-1))
               ff=vol(I)/vol(I-1)
               if (ff.gt.1) ff=1/ff

            endif
            ff=1
            if (I.lt.NS) then
               dh=dabs( phshell(I)-phshell(I+1))
               dh1=0d0
               if (i.ge.2) dh1=dabs( phshell(I)-phshell(I-1))
               if (dh1.gt. dh) dh= dh1
               daw=dabs( awshell(I)-awshell(I+1))
               ff=vol(I)/vol(I+1)
               if (ff.gt.1) ff=1/ff
        endif

        if (ff.le. 0.5d0 .and.daw.le..02 .and.dh.le..04  ) ismerge=1

        xnca=0d0
             if(ismerge.eq.1 ) then
                npl=np-npsolid
               
                DO II=1,NS
         xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
                        print*, 'before merge ', II, phshell(II)

                     enddo
                     print*,'befor merge ca', xnca
                     
             if (I.eq.NS) then
                  DO J=1,1,NP
                     xn(NS-1,j)=xn(NS-1,j)+xn(NS,j)
                  enddo

                  DO J=1,1,NP
                     xn(NS,j)=0d0
                  enddo

                  NS=NS-1
                  x(NS+1)=x(NS+2)
                print*, ' The shell'  , I, 'is merged to, shell',  I-1

                  
                  else

             print*, ' The shell'  , I, 'is merged to, shell',  I+1,ns
                     DO J=1,NP
                     dd=xn(I,j)+xn(I+1,j)
                     xn(I,j)=dd
                  enddo
                  NS=NS-1

                  x(I+1)= x(I+2)
                  
                  DO Ii=I+1,NS
                         DO J=1,NP
                            xn(Ii,j)=xn(Ii+1,j)
                         enddo
                     x(Ii+1)= x(Ii+2)
                      enddo
                      
                  DO J=1,NP
                            xn(NS+1,j)=0
                  enddo
                  
                  
              
                  
                     endif

                     xnca=0
                     npl=np-npsolid
                     
             DO II=1,NS
         xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         

                ml(1)=1000/mm(1)
                
               DO J=2,np
                  ml(j)= ml(1)*xn(II,J)/xn(II,1)

               enddo
               IS=II
               DO kk=1,10
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           enddo
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
          awshell(II)= aw
               
           write(6,'(A,I5,5E15.6)')'After merge ', Ii,aw,phshell(Ii),ta
            enddo
            print*,'after merge ca' ,xnca
            return
                   endif

             
          enddo


          return
          end

c     Redistribuete starting from shell 1

      subroutine merge_redis(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np),xna(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),ml(np),flux(NSMM,np+npsolid)
      
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell
     & ,vshell0(nsmm)
c     common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa

      real*8 dx1(nsmm)
      
       common /dxmin/dxmin0
      common /NS0/NS0

      if (NS.le.1) return
       if(ns.gt.ns0) return
       
      vmax=0
      DO I=1,NS
         vol(i)= xn(I,1)*mv(1)+xn(I,2)*mv(2)
         DO J=6,np-npsolid
         vol(I)=vol(I)+xn(I,j)*mv(j)
         enddo
      enddo

      vmax=0
      
      dxmin=dxmin0*.75
c      print*, dxmin, dxmin0
      NS1=ns
      DO JJ=1,ns1
c      if (NS.le.ns0) return

      rcore3= x(1)**3
      DO I=ns,1,-1
         
         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( rcore3+vs/4/pi*3)**(1/3d0)

         vdry= vs+ xn(I,2)*mv(2)+mv(1)*xn(I,1)

       DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rr= ( rcore3+vdry/4/pi*3)**(1/3d0)
         dx=rr-xss
         dx1(I)=dx
      enddo
      
c         print*,' dx, dxmin', I,dx, dxmin
      
      DO I =1,ns



         

         if (dx1(I).le.dxmin) then

            I1= I


            if (I.gt.1.and.I.lt.NS) then
               if (vol(I-1).ge.vol(I+1)) I1=I-1
            endif
            if (I1.ge.ns) I1=Ns-1
            
            
            if(NS.le.ns0) then
           print*, 'merge  redistribute ', I1,dx
           xna0=0
           DO II=1,NS
              xna0=xna0+xn(II,16)
           enddo
           DO J=1,np-npsolid
              ss=xn(I1,j)+xn(I1+1,j)
              xn(I1,j)=ss/2
              xn(I1+1,j)=ss/2
           enddo
           write(6,'(A,100E15.6)')'merge redis inner',(dx1(kk),kk=1,ns)
           
            goto 222
               
            endif
            
            return
            
           print*, 'merge dx ', ns,dx
           
            
            DO J=1,np
               xn(I1,J)=xn(I1,J)+xn(I1+1,J)
            enddo
            X(i1+1)= X(i1+2)

            ns=ns-1
            do iI=i1+1,NS
               do j=1,nP
               xn(II,J)=xn(II+1,J)
            ENDDO
            X(ii+1)= X(ii+2)
            ENDDO
            
            goto 222
         endif
         
c     5*dx
         dd=x(NS+1)-x(1)/10
c         print*,'ns, dd,dx,', ns,dd,dx
         
         if (NS.ge.10) then
            if (dx.le. dd) then
             if (I.gt.1) then

                dh=dabs( phshell(I)-phshell(I-1))
               if (i.lt.NS) dh1=dabs( phshell(I)-phshell(I+1))
               dh1=0
               if (dh1.gt. dh) dh= dh1
            else
              dh=dabs( phshell(I)-phshell(I+1))
               
           endif
           
             if (I.gt.1) then
                daw=dabs( awshell(I)-awshell(I-1))
                daw1=0d0
                if (i.lt.NS) daw1=dabs( awshell(I)-awshell(I+1))

               if (daw1.gt. daw) daw= daw1
            else
              daw=dabs( awshell(I)-awshell(I+1))
               
            endif

            
c
            print*,'dh, daw ', I,dh,daw
            
            if (dh.le. 0.05 .and.daw.le..01) then
               print*, 'merge, 5Xdxmin0m', I, dx, dh
            I1= I
            if (I.eq.NS) I1=I-1
            if (I.eq.1) I1=1
            if (I.gt.1.and.I.lt.NS) then
               I1= I
               if (vol(I-1).lt.vol(I+1)) I1=I-1
            endif
            

            DO J=1,np
               xn(I1,J)=xn(I1,J)+xn(I1+1,J)
            enddo
            X(i1+1)= X(i1+2)

            ns=ns-1
            do iI=i1+1,NS
               do j=1,nP
               xn(II,J)=xn(II+1,J)
            ENDDO
            X(ii+1)= X(ii+2)
            ENDDO
            
            goto 222


            endif
            

            endif
         endif
         

      enddo
      return
 222  continue
                     xnca=0
             DO II=1,NS
         xnca=xnca+ xn(Ii,16)
      enddo
      
         npl=np-npsolid
c                     timea=1
c                     call cal_flux(timea,x,flux,dtime,NS)
                     
             DO II=1,NS


                ml(1)=1000/mm(1)
                
               DO J=2,np
                  ml(j)= ml(1)*xn(II,J)/xn(II,1)

               enddo
c               DO KK=1,10
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

        
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
          awshell(II)= aw
c          print*,ii, aw,x(II+1)-x(II)
       enddo

c            enddo
            print*,'after merge d x na' ,xnca,xna0
            return

         enddo

      if (NS.le.ns0) return
         
         DO I=1,NS
c     compare shell below 
      ismerge =0

             if (I.gt.1) then
               dh=dabs( phshell(I)-phshell(I-1))
               if (i.lt.NS) dh1=dabs( phshell(I)-phshell(I+1))
               if (dh1.gt. dh) dh= dh1
               daw=dabs( awshell(I)-awshell(I-1))
               ff=vol(I)/vol(I-1)
               if (ff.gt.1) ff=1/ff

            endif
            ff=1
            if (I.lt.NS) then
               dh=dabs( phshell(I)-phshell(I+1))
               dh1=0d0
               if (i.ge.2) dh1=dabs( phshell(I)-phshell(I-1))
               if (dh1.gt. dh) dh= dh1
               daw=dabs( awshell(I)-awshell(I+1))
               ff=vol(I)/vol(I+1)
               if (ff.gt.1) ff=1/ff
        endif

        if (ff.le. 0.5d0 .and.daw.le..02 .and.dh.le..04  ) ismerge=1

        xnca=0d0
             if(ismerge.eq.1 ) then
                npl=np-npsolid
               
                DO II=1,NS
         xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
                        print*, 'before merge ', II, phshell(II)

                     enddo
                     print*,'befor merge ca', xnca
                     
             if (I.eq.NS) then
                  DO J=1,1,NP
                     xn(NS-1,j)=xn(NS-1,j)+xn(NS,j)
                  enddo

                  DO J=1,1,NP
                     xn(NS,j)=0d0
                  enddo

                  NS=NS-1
                  x(NS+1)=x(NS+2)
                print*, ' The shell'  , I, 'is merged to, shell',  I-1

                  
                  else

             print*, ' The shell'  , I, 'is merged to, shell',  I+1,ns
                     DO J=1,NP
                     dd=xn(I,j)+xn(I+1,j)
                     xn(I,j)=dd
                  enddo
                  NS=NS-1

                  x(I+1)= x(I+2)
                  
                  DO Ii=I+1,NS
                         DO J=1,NP
                            xn(Ii,j)=xn(Ii+1,j)
                         enddo
                     x(Ii+1)= x(Ii+2)
                      enddo
                      
                  DO J=1,NP
                            xn(NS+1,j)=0
                  enddo
                  
                  
              
                  
                     endif

                     xnca=0
                     npl=np-npsolid
                     
             DO II=1,NS
         xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         

                ml(1)=1000/mm(1)
                
               DO J=2,np
                  ml(j)= ml(1)*xn(II,J)/xn(II,1)

               enddo
               IS=II
               DO kk=1,10
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           enddo
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
          awshell(II)= aw
               
           write(6,'(A,I5,5E15.6)')'After merge ', Ii,aw,phshell(Ii),ta
            enddo
            print*,'after merge ca' ,xnca
            return
                   endif

             
          enddo


          return
          end


      
      subroutine split_shells(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP),ml(np) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),xs(NSMM+1),dxliqdry(NSMM)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),vshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell,vshell
c       common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

       common /dxmin/dxmin0
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa

      vcore= 4*pi/3*x(1)**3
      rcore3= x(1)**3
      DO I=1,ns

         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( rcore3+vs/4/pi*3)**(1/3d0)
         vdry= vs+ xn(I,2)*mv(2)
         DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rdry= (rcore3+ vdry/4/pi*3)**(1/3d0)
         dxliqdry(I)=rdry-xss
      enddo
      dxmin=dxmin0*1.5
      
      

c     not possible to merge
      
      if (NS.ge.50) return
      I1=0
      if (NS.eq.2) then
         I=NS

         dx= dxliqdry(I-1)
            if (dxliqdry(I). Lt. dxliqdry(I-1) ) dx= dxliqdry(I)

         dhmax=.5
         if (dx.le.dxmin*3) dhmax=1d0

         dh=dabs(phshell(I)-phshell(I-1))
         if (dh.ge.dhmax) then
            I1=I-1
            if (dxliqdry(I). gt. dxliqdry(I-1) ) I1=I
         endif
            dx= dxliqdry(I1)
            if (dx.le.dxmin) I1=0

            if (dx.le.1.5*dxmin0 .and. I1.eq.2) I1=0
            
            if (I1.ge.1) goto 333
         
      endif

      
      DO I=NS-1,2,-1
c         print*,'s', I, phshell(I),phshell(I+1)
         dh=dabs(phshell(I)-phshell(I-1))
         dx= dxliqdry(I-1)
            if (dxliqdry(I). Lt. dxliqdry(I-1) ) dx= dxliqdry(I)

         dhmax=.5
         if (dx.le.dxmin*3) dhmax=1d0

         if (dh.ge.dhmax) then

            I1=I-1
            if (vshell(I).gt. vshell(I-1)) I1=I
            dx= dxliqdry(I1)
            if (dx.le.dxmin) I1=0
            if (I1.eq.ns .and. dx.le.2*dxmin) I1=0
            
            if (I1.ge.1) goto 333
            
            
         endif

         I1=0
         
         dx= dxliqdry(I+1)
            if (dxliqdry(I). Lt. dxliqdry(I+1) ) dx= dxliqdry(I)

         dhmax=.5
         if (dx.le.dxmin*3) dhmax=1d0
         
            dh=dabs(phshell(I)-phshell(I+1))
            
            if (dh.ge.dhmax)  then
            I1=I+1
            if (vshell(I).gt. vshell(I+1)) I1=I

            endif

            dx= dxliqdry(I1)
            if (dx.le.dxmin) I1=0
c     fort outermost shell no split
            if (I1.eq.ns .and. dx.le.2*dxmin) I1=0
            
            if(I1.gt.0) then
         goto 333
         endif

      enddo
 333  continue

      if (I1.ge.1) then
         print*, 'split I', I1
         xnca=0
            npl=np-npsolid
            DO II=1,NS
c               print*,'before split ', Ii,phshell(Ii)
         xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
            enddo
            print*, 'before ca', xnca
         DO I2=NS,I1+1,-1
            x(I2+2)=x(I2+1)
            DO J=1,Np
               xn(I2+1,j)= xn(I2,j)
            enddo
         enddo

         DO J=1,NP-npsolid
            dd=xn(I1,j)/2
            xn(I1,j)= dd
            xn(I1+1,j)= dd
         enddo

          DO J=NP-npsolid+1,np
            xn(I1+1,j)= 0d0
         enddo

         dx= x(I1+1)-x(I1)
         x(I1+2)= x(I1+1)
         x(I1+1)= x(I1)+dx/2
         NS=NS+1
         xnca=0
            DO II=1,NS
               ml(1)=1000/mm(1)
               DO J=2,np
                  ml(J)= ml(1)/xn(ii,1)*xn(ii,j)
               enddo
               IS=Ii
               DO kk=1,10
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
        enddo
        
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
          awshell(II)= aw
               
c           print*,'after split ', Ii,phshell(Ii)
           xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
          enddo
          print*,'after ca', xnca
          
      
       endif
       

          return
          end

c     split aw
            subroutine split_shells_aw(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP),ml(np) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),xs(NSMM+1),dxliqdry(NSMM)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),vshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell,vshell

c       common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

       common /dxmin/dxmin0
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa
 
      
      DO I=1,ns

         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( vs/4/pi*3)**(1/3d0)
         vdry= vs+ xn(I,2)*mv(2)
         DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rdry= ( vdry/4/pi*3)**(1/3d0)
         dxliqdry(I)=rdry-xss
      enddo

      dxmin=dxmin0*1.5
      
c     not possible to merge
      
      if (NS.ge.50) return
      I1=0
      DO I=NS-1,2,-1
c         print*,'s', I, phshell(I),phshell(I+1)
         dh=dabs(awshell(I)-awshell(I-1))
         if (dh.ge..1) then

            I1=I-1
            if (vshell(I).gt. vshell(I-1)) I1=I
             dx=dxliqdry(I)

            if (dx.le.dxmin) I1=0
            if (I1.eq.ns .and. dx.le.2*dxmin) I1=0
            
         endif

         
            dh=dabs(awshell(I)-awshell(I+1))
            
            if (dh.ge..1)  then
            I1=I+1
            if (vshell(I).gt. vshell(I+1)) I1=I


            endif

            dx= dxliqdry(I1)
            if (dx.le.dxmin) I1=0
            if (I1.eq.ns .and. dx.le.2*dxmin) I1=0
            
            if(I1.gt.0) then
            write(6,'(A,I5,5E15.6)') 'split I aw', I1,
     *awshell(I1-1), awshell(i1),awshell(i1+1)
            xnca=0
            npl=np-npsolid
            DO II=1,NS
               print*,'before split aw', Ii,awshell(Ii)
           xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
          enddo
          print*,'bfor ca', xnca

         DO I2=NS,I1+1,-1
            x(I2+2)=x(I2+1)
c     mve I1+1 to I1+2
            DO J=1,Np
               xn(I2+1,j)= xn(I2,j)
            enddo
         enddo
         
         DO J=1,NP-npsolid
            xn(I1,j)= xn(I1,j)/2
            xn(I1+1,j)= xn(I1,j)
         enddo
         DO J=NP-npsolid+1,np
            xn(I1+1,j)= 0d0
         enddo

         dx= x(I1+1)-x(I1)
         x(I1+2)= x(I1+1)
         x(I1+1)= x(I1)+dx/2
         NS=NS+1
         xnca=0d0
            DO II=1,NS
               ml(1)=1000/mm(1)
               DO J=2,np
                  ml(J)= ml(1)/xn(ii,1)*xn(ii,j)
               enddo
               IS=ii
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         
          awshell(II)= aw
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
               
           print*,'after split aw', Ii,awshell(Ii)
           xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
          enddo
          print*,'after aw ca', xnca
            

         return
         endif

      enddo
      


          return
          end


      subroutine split_shells_ReSAM(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP),ml(np) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np)
      real*8 x(NSMM+1),dxliqdry(NSMM)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),vshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell,vshell
c       common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

       common /dxmin/dxmin0
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa

      
      DO I=1,ns

         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( vs/4/pi*3)**(1/3d0)
         vdry= vs+ xn(I,2)*mv(2)
         DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rdry= ( vdry/4/pi*3)**(1/3d0)
         dxliqdry(I)=rdry-xss
      enddo
      dxmin=dxmin0*1.5
      
      

c     not possible to merge
      
      if (NS.ge.50) return
      I1=0
      DO I=NS-1,2,-1
c         print*,'s', I, phshell(I),phshell(I+1)
         dh=dabs(phshell(I)-phshell(I-1))
         if (dh.ge..5) then

            I1=I-1
            if (vshell(I).gt. vshell(I-1)) I1=I
            
         endif

         
            dh=dabs(phshell(I)-phshell(I+1))
            
            if (dh.ge..5)  then
            I1=I+1
            if (vshell(I).gt. vshell(I+1)) I1=I


            endif

            dx= dxliqdry(I1)
            if (dx.le.dxmin) I1=0
            
            if(I1.gt.0) then
            print*, 'split I', I1
      xnca=0
            npl=np-npsolid
            DO II=1,NS
               print*,'before split ', Ii,phshell(Ii)
         xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
            enddo
            print*, 'before ca', xnca
         DO I2=NS,I1+1,-1
            x(I2+2)=x(I2+1)
            DO J=1,Np
               xn(I2+1,j)= xn(I2,j)
            enddo
         enddo

         DO J=1,NP-npsolid
            dd=xn(I1,j)/2
            xn(I1,j)= dd
            xn(I1+1,j)= dd
         enddo

          DO J=NP-npsolid+1,np
            xn(I1+1,j)= 0d0
         enddo

         dx= x(I1+1)-x(I1)
         x(I1+2)= x(I1+1)
         x(I1+1)= x(I1)+dx/2
         NS=NS+1
         xnca=0
            DO II=1,NS
               ml(1)=1000/mm(1)
               DO J=2,np
                  ml(J)= ml(1)/xn(ii,1)*xn(ii,j)
               enddo
               is=II
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
          awshell(II)= aw
               
           print*,'after split ', Ii,phshell(Ii)
           xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
          enddo
          print*,'after ca', xnca
          
         return
         endif

      enddo
      


          return
          end

c     split aw
            subroutine split_shells_aw_resam(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP),ml(np) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np)
      real*8 x(NSMM+1),dxliqdry(NSMM)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),vshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell,vshell
c       common/flux/T,Ta
       common /dxmin/dxmin0
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac


      
      DO I=1,ns

         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( vs/4/pi*3)**(1/3d0)
         vdry= vs+ xn(I,2)*mv(2)
         DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rdry= ( vdry/4/pi*3)**(1/3d0)
         dxliqdry(I)=rdry-xss
      enddo

      dxmin=dxmin0*1.5
      
c     not possible to merge
      
      if (NS.ge.50) return
      I1=0
      DO I=NS-1,2,-1
c         print*,'s', I, phshell(I),phshell(I+1)
         dh=dabs(awshell(I)-awshell(I-1))
         if (dh.ge..1) then

            I1=I-1
            if (vshell(I).gt. vshell(I-1)) I1=I
            
         endif

         
            dh=dabs(awshell(I)-awshell(I+1))
            
            if (dh.ge..1)  then
            I1=I+1
            if (vshell(I).gt. vshell(I+1)) I1=I


            endif

            dx= dxliqdry(I1)
            if (dx.le.dxmin) I1=0
            
            if(I1.gt.0) then
            write(6,'(A,I5,5E15.6)') 'split I aw', I1,
     *awshell(I1-1), awshell(i1),awshell(i1+1)
            xnca=0
            npl=np-npsolid
            DO II=1,NS
               print*,'before split aw', Ii,awshell(Ii)
           xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
          enddo
          print*,'bfor ca', xnca

         DO I2=NS,I1+1,-1
            x(I2+2)=x(I2+1)
c     mve I1+1 to I1+2
            DO J=1,Np
               xn(I2+1,j)= xn(I2,j)
            enddo
         enddo
         
         DO J=1,NP-npsolid
            xn(I1,j)= xn(I1,j)/2
            xn(I1+1,j)= xn(I1,j)
         enddo
         DO J=NP-npsolid+1,np
            xn(I1+1,j)= 0d0
         enddo

         dx= x(I1+1)-x(I1)
         x(I1+2)= x(I1+1)
         x(I1+1)= x(I1)+dx/2
         NS=NS+1
         xnca=0d0
            DO II=1,NS
               ml(1)=1000/mm(1)
               DO J=2,np
                  ml(J)= ml(1)/xn(ii,1)*xn(ii,j)
               enddo
           call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         
          awshell(II)= aw
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
               
           print*,'after split aw', Ii,awshell(Ii)
           xnca=xnca+ xn(Ii,36)+ xn(Ii,npl+10)+xn(Ii,npl+11)         
          enddo
          print*,'after aw ca', xnca
            

         return
         endif

      enddo
      


          return
          end



c     --------------------------------------------------------
c
c     Calculate the Equilibrium composition of solution with NaCl crystal
c     Ma: the inital composition as input and output is the composition with NaCl Crystal
c     Msalt: is the molalty of NaCl crystal
c      
c     -------------------------------------------------------

      
      subroutine effl_misch(Ta,ma,msalt)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50)
        real*8 fcnmisch
	external  fcnmisch
        real*8 MA(NP),m(NP),msalt
	common/misch/ T, M
        common /kout/xx2

        
        
        t=ta
              do kk=1,np

        M(kk)= MA(kk)
        enddo
 	xmin=1D-14              
          	 xmax=m(16)-1D-20        
         xmax1=m(17)-1D-20        
        if (xmax1.le.xmax) xmax= xmax1
        
	erabs=0.d0
	errel=0d0
	ITmax=2000

	xx1= fcnmisch(xmin)
	xx2= fcnmisch(xmax)


        



	if(xx1.le.0.) then
	xmax=xmin
	goto 2
      endif
      
	if(xx2.ge.0.) 	goto 2

c       print*, 'xmin effmish', xmin,xx1
c       print*, 'xmax effmisch', xmax,xx2
	call dzbrens(fcnmisch, erabs,errel,xmin,xmax,ITMAX)
 2      continue

        Msalt= xmax
        MA(16)=Ma(16)-xmax
        MA(17)=Ma(17)-xmax

      call aw_back(t,Ma,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
              Ap= Ma(16)*ma(17)*gammaCl*gammaNa


              print*, ap
              


        
 999            return
	end  

c     ----------------------------------
c    Required by eff_misch
c
c     ---------------------------------
           function fcnmisch(x)

        IMPLICIT REAL*8 (A-H,m,O-Z)
      integer NP
      parameter (np=50)
      real*8 M(NP), NL(NP),msalt,ML(NP)
      common/misch/ T, M

      
      
      msalt= x
      do kk=1,np

         ML(kk)=M(kk)
         enddo
      
      ML(16)=ML(16)-Msalt
      if (ML(16).le.1D-20) ML(16)=1D-20
      
      ml(17)=ML(17)-Msalt
      if (ML(17).le.1D-20) ML(17)=1D-20



      
      call aw_back(t,Ml,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNA)

       Aproduct= ML(16)*ml(17)*gammaCl*gammaNa
      
       fcnmisch= Aproduct -Apnacl(T)
       

      
         
        return
	end


c     apNH42C2O4 in hPa**3 of NH42_C2O4, obtained from the solubitlity data
c     Hefter 2018


      function apnahc2o4(T)
        IMPLICIT REAL*8 (A-H,m,O-Z)


        apnahc2o4 = 8D-6
      return
      end

      function apna2c2o4(T)
        IMPLICIT REAL*8 (A-H,m,O-Z)
        apna2c2o4 =  0.21329E-01 *dexp (-0.17274E+04*(1/t-1/298.15))


      return
      end



c     Hill 1935
      function apNH42C2O4(T) 
        IMPLICIT REAL*8 (A-H,m,O-Z)
       	fcal=-0.39738E+02
     &  -0.24393E+05*(1/t-1/298.15d0)+ 0.49090E+02* dlog(t/298.15d0)
        apNH42C2O4= dexp(fcal)
      return
      end


c      vapour pressure product NH4HC2O4 + 05H2O   in hPa**2
c otained from solubilty data,  Buttke LG et al 2016 
      function apNH4HC2O4_05h2o(T)
        IMPLICIT REAL*8 (A-H,m,O-Z)
       	fcal=  -0.30716E+02  -0.24169E+05*(1/t-1/298.15d0)
        apNH4HC2O4_05h2o= dexp(fcal)
      return
      end




c     --------------------------------------------------
c
c     calculate the flux f1 (from the shell to to crystal) 
c      and f2 (the flux liquid to crystal of the same shell)
c
c     NS: the total number of shells
c      II: the location of the crystal 
c     
c     --------------------------------------------------
c     ammonium oxalate crystal 
c 

      subroutine  cal_ao_flux(NS, iI,Tdrop, x,f1,f2)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      parameter (np=50,NSMM=100,npsolid=12)
      real*8  x(*),ml(NP),xa(NSMM),ml0(np),ml1(NSMM),mlm1(NSMM)
      
       
      real*8 dl_factor2(2,np),dl_factor(NP)
     
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin

      real*8 xn(NSMM,NP),x1(NP)
      common /xn/xn
      common/eutectic/ feutectic0


      common /output/imode_output,idiff,imode_pH,imode_eq
      
      common /ienhance/Ienh,iscenter,isAHO

      common /pi/pi
      
      
       common/flux/Taa,Taaa,pressaaa,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &      partvaphcl, partvapOA,partvaplac

c        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
c     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak



      real*8 MM(NP) ! molar Mass
      real*8 Mv(NP) ! molar volume
      
      common /m/ mm, mv,izc(NP)
      common /Nacl/T,ML,ML0

      common /testA/ Ierr

       common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseq

       
       real*8 gamma2(NSMM*2,NP)
      common /gamma2/gamma2
      integer NSP(npsolid)
      integer NSP_index(npsolid,10)
      real*8 xNSP_nv(npsolid,10)
      COMMON /nsp/ NSP, NSP_INDEX
      COMMON /xnsp/ xNSP_NV
      common /xvion/ xvion
      
      feutectic=feutectic0
      

      

      

       fac=.5


      Ierr=0

      I=II

      f1=0d0

      xvsolid= 0d0
      DO kk=1,npsolid
         iss= np-npsolid+kk
      xvsolid=xvsolid+mv(iss) *xn(I,iss)
      enddo

      vol= 4*pi/3d0*x(I)**3 +xvsolid

      xvliquid =4*pi/3d0*( x(I+1)**3-x(I)**3) -xvsolid
      r1=(vol/4d0/pi*3)**(1/3d0)
      r2=x(I+1)

      if (I.lt.NS) then
         if (r2.le. x(I+2)/2) r2 = x(I+2)/2
      endif
      
      v1= xvliquid
      
            
      ML(1) = 1000d0/mm(1)

      t=Tdrop

c     take the mean composition of shell I and I+1
      DO J=1,NP
         ML1(J)=ML(1)* xn(I,J)/xn(I,1)
                  if (I.lt.NS) then
         ML(j)= 
     *(ML(1)* xn(I,J)/xn(I,1)*fac+(1-fac)*ML(1)* xn(I+1,J)/xn(I+1,1))
                 else
         ML(j)= ML(1)* xn(I,J)/xn(I,1)
         ML1(J)=ML(1)* xn(I,J)/xn(I,1)
      endif
c         print*, i,j, ML(j) 
      enddo

      do kk=1,np
      ML0(kk)=ML(kk)
      enddo
c     NH4HC2O4 .05H2O
      is=I      
      IP=ISAHO
      spp1=ss_solid(isAHO, T,ml1)
            spp=ss_solid(isAHO, T,ml)
            xnv=0d0
            DO K=1, nsp(IP)
               j= nsp_index(Ip,k)
               if (j.ge.2) xnv=xnv+ xnsp_nv(IP,k)
            enddo

            spp=spp !**(1d0/xnv)             !/xnv
            
      Ienh=0
      if (iscenter.eq.1) ienh=1

      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        



      
      do kk=1,np
       ML0(kk)=ML(kk)
       enddo

          ff=spp
          if (ff.le.0.01) ff=0.01

          


      
      
c     Na2C2O4
      fs=0
c     it is controlled by the species with smallest concentration
      cmin=1D10

      
      DO K=1,NSP(IP)
               j= nsp_index(Ip,k)
               if(J.eq.1) then

                  call caldl(t,aw,ml,x,J,dl)
               endif

               if(J.gt.1) then
       c1 = xn(I,1)/v1* ML(j)/ML(1)
       if (c1.lt.cmin) then 
          cmin=c1
c     newdl
         call caldl(t,aw,ml,x,J,dl)

          

          dln= dl* dl_factor2(1,j)
       endif
      endif
      
      enddo
      c1=cmin
      c0 =cmin/ff
       f2 = dlN  * (c1-c0)*4 *pi*r1*(r2)/dabs(r2-r1)
      fspp=ff
c      write(6,'(A, 2I5,6E15.6)') 'FF', I, IP, dln
c     &     , c1,c0, r1,r2, dabs(r2-r1)
      
      
c     consider the eutectic structure
          if (feutectic.gt.0d0) then
             dd=dabs(r2-r1)
c     consider the gradient
             ff=dd/feutectic
             f2=f2*ff
c     number of layers
             ff=dd/feutectic
c                         if (ff.le.1) ff=1, if layer is 
            f2=f2*ff
         endif
          

         if ((Ip.eq.10 .or. ip.eq.11) .and.I.eq.1) then
c            write(6,'(A,2I5,5E15.6)')'ff2', I,Ip, f2,spp,fspp
         endif

         
c        write(39,'(A,11E15.6)') 'spp ', spp1, spp, f2

c          print*, 'f2 ', f2


       
       if (spp1.le.1d0 .and. f2.gt.0d0) f2=0d0 ! no growth when spp1 <1
       if (spp1.gt.1d0 .and. f2.lt.0d0) f2=0d0 ! no decrease when spp1 >1



       
       f1=0d0
       
       if (I.ge.2) then
          
       v2 = 
     & 4*pi/3 * x(I-1)**3 

       DO  kk=1,npsolid
       iss=np-npsolid+kk
       v2=v2+mv(iss)*xn(I-1,iss) 
       enddo

      r2=(v2/4/pi*3)**(1/3d0)  ! inner radius of shell I-1
      r1=x(I)  ! outer radius I-1
      if (r2.le.1D-30) r2= r1/2
      
      
      vliq= xn(I-1,1) *mv(1)+ xn(I-1,2)*mv(2)
      DO J= 6,NP-npsolid
         vliq= vliq + xn(I-1,J)*mv(J)
      enddo
      v1=vliq

         

c     take the mean composition of the shell I and I+1
      DO J=1,NP
         I1= I-1
         ML(j)= ML(1)* xn(I1,J)/xn(I1,1)
         ML1(j)= ML(1)* xn(I1,J)/xn(I1,1)
         if (i.ge.3)
     * ML(j)= 
     *(ML(1)* xn(I1-1,J)/xn(I1-1,1)*(1-fac)+ 
     *fac*ML(1)* xn(I1,J)/xn(I1,1))
      enddo

      IS=I1
c     NH4HC2O4 .05H2O
      call vapnew(T,ML1,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
      spp1=ss_solid(isAHO, T,ml1)
            spp=ss_solid(isAHO, T,ml)
c            spp=spp**(1 /xnv)
                 ff=spp
           if (ff.le.0.01d0) ff=0.01
      cmin=1D10
      
      DO K=1,NSP(IP)
               j= nsp_index(Ip,k)
               if(J.gt.1) then
         call caldl(t,aw,ml,x,J,dl)
                  endif
                  if(J.gt.1) then
       c1 = xn(I,1)/v1* ML(j)/ML(1)
       if (c1.lt.cmin) then

         call caldl(t,aw,ml,x,J,dl)


          cmin=c1
                  dln= dl* dl_factor2(1,j)
       endif
      endif
      
      enddo
      c1=cmin
      c0 =cmin/ff
       f1 = dlN  * (c1-c0)*4 *pi*r1*(r2)/dabs(r2-r1)



c     consider the eutectic structure
          if (feutectic.gt.0d0) then
             dd=dabs(r2-r1)
c     consider the gradient
                         ff=dd/feutectic
                         if (ff.le.1) ff=1
            f1=f1*ff
         endif
              

       if (spp1.le.1d0 .and. f1.gt.0d0) f1=0d0 ! no growth when spp1 <1
       if (spp1.gt.1d0 .and. f1.lt.0d0) f1=0d0 ! no decrease when spp1 >1



       endif

      
       
        return
        end



c     determine the equilibroum compostion of NH42 oxalate



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NEW this is switched off!!!!!!!!!!!!!!!!
c      
c     --------------------------------------------------
c
c     Required by cal_AO_flux, to calculate the equilibrium NH4+ and OA-- concentration over ammonium oxalate Crystal
c      
c     --------------------------------------------------
c     doen't work for Na ML(16)









c     --------------------------------------------------------
c
c     Calculate the Equilibrium composition of solution with NaCl crystal
c     Ma: the inital compostion as input and output is the composition with NaCl Crystal
c     Msalt: is the molalty of oxalate crystal
c      
c     -------------------------------------------------------


c     NH4HC2O4 05H2O
      subroutine effl_AO(Ta,ma,msalt)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50,NSMM=100,npsolid=12)
        real*8 fcnmisch
	external  fcnsolidao
        real*8 MA(NP),m(NP),msalt
	common/efflAO/ T, M
        common /kout/xx2
        common /test/ Ierr
       common /Ienhance/Ienh,iscenter,isAHO

       
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

       real*8 gamma2(NSMM*2,NP)
       common /gamma2/gamma2

       
c        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
c     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak

        integer NSP(npsolid)
      integer NSP_index(npsolid,10)
      real*8 xNSP_nv(npsolid,10)
      COMMON /nsp/ NSP, NSP_INDEX
      COMMON /xnsp/ xNSP_NV

        
        Ierr=0
        
        
        t=ta
        DO kk=1,np
        M(kk)= MA(kk)
        enddo
        msalt=0d0
        ss= ss_solid(ISAHO, t, ma)

        if (ss.le.1d0) return
           
 	xmin=0d0
        xmax=100d0
        DO k=1, nsp(ISAHO)
           j= nsp_index(Isaho,k)
           if (j.eq.14) then
              ma(14) =ma(14)+ma(15)
              ma(15)=0d0
              m(14)=ma(14)
              m(15)=ma(15)
           endif
           
           if (J.ge.2 .and. ma(j).le.xmax) xmax= ma(J)
        enddo
        
c        print*,'xmin', xmin, xmax
c        print*, m(16),m(26)
        
	erabs=0.d0
	errel=0d0
	ITmax=200

	xx1= fcnsolidAO(xmin)
	xx2= fcnsolidAO(xmax)
c
	write(6,'(A, 4E14.6)')'effl_ao', xmin,xmax,xx1,xx2
	if(xx1.le.0.) then
	xmax=xmin
	goto 2
      endif
      
	if(xx2.ge.0.) 	goto 2

c       print*, 'xmin ao', xmin,xx1
c       print*, 'xmax ao', xmax,xx2
	call dzbrens(fcnsolidAO, erabs,errel,xmin,xmax,ITMAX)
 2      continue
	xx2= fcnsolidAO(xmax)
c        print*, 'xx2 ', xx2

        if (ierr.ge.1) then
           print*, 'error: eff_AHO'
           stop

           endif


        Msalt= xmax
        DO k=1, nsp(ISAHO)
           j= nsp_index(Isaho,k)
           if(j.ge.2)    ma(J)=ma(j)-xnsp_nv(isaho,k)*msalt
        enddo

      

      
      call vapnew(T,Ma,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)




 999            return
	end  




c     isAHO=1 : NH4HC2O4
c     else      (NH4)2C2O4

           function fcnsolidAO(x)

        IMPLICIT REAL*8 (A-H,m,O-Z)
      integer NP
      parameter (np=50,NSMM=100,npsolid=12)
      real*8 M(NP), NL(NP),msalt,ML(NP)
      common/efflAO/ T, M
      common /Ienhance/Ienh,iscenter,isAHO

       
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

       real*8 gamma2(NSMM*2,NP)
       common /gamma2/gamma2

      
c        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
c     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak

      integer NSP(npsolid)
      integer NSP_index(npsolid,10)
      real*8 xNSP_nv(npsolid,10)
      COMMON /nsp/ NSP, NSP_INDEX
      COMMON /xnsp/ xNSP_NV
      
      
      msalt= x
c      ML=M
      do kk=1,np
         ml(kk)=m(kk)
      enddo
      

      DO k=1, nsp(ISAHO)
           j= nsp_index(Isaho,k)
           if(j.ge.2)    ml(J)=ml(j)-xnsp_nv(isaho,k)*msalt
        enddo

        ss= ss_solid(isaho,t,ml)
        
        fcnsolidAO=ss-1

        return
	end



      

c     -----------------------------------------------------------------------------

c

      	subroutine calHNew(Ta,ma)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50,NSMM=100, npsolid=12)
        real*8 fcnnew
	external  fcnnew

        
        real*8 MA(NP),m(NP),xn(NSMM, np)
        common /xn/xn

	common/suls/ T, M, ppartco2

       
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

       common /kout/xx2
       
      common /output/imode_output,idiff,imode_pH,imode_eq

       common /time/ time, xn_area
       common /test/ Ierr
       real*8 gamma2(NSMM*2,NP)
      common /gamma2/gamma2
c      print*,'is in calhnew',is
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     &     ,gammaH2OA       ,gammapo4
        
        
      if (IS.gt.2*NS) IS =1
      if (ns.lt.1) IS=1
      if (Is.lt.1) IS=1

      ii=0
      DO j=1,NP
c         gamma2(Is,j)=1d0
      enddo
      gammah0 = gamma2(IS,6)

      
       DO iiI=1,10
       
c interat H+
         ii=0

         call calHNewaa(Ta,ma)

         call aw_back(ta,Ma,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         I=IS
       gamma2(I,13) =get_gammaco2(ta,Ma)       
       gamma2(I,6)=gammah
       gamma2(I,12)=gammaNH4
       gamma2(I,17)=gammaCl
       gamma2(I,16)=gammaNa
       gamma2(I,18)=gammaNO3
       gamma2(I,27)=gammas1
       gamma2(I,28)=gammaS2

       gamma2(I,8)=gammaNO3
       gamma2(I,34)=gammaNO3
       
       gamma2(I,1)=aw
       xmplus=ma(6)+ma(19)+ma(12)+ma(16)
       gammaOH =get_gammaOH(xMplus)
       gamma2(I,7)=gammaoH

       gamma2(I,15)=gammaHCO3
       gamma2(I,14)=gammaCO3

       gamma2(I,24)=gammaH3PO4
       gamma2(I,25)=gammaH2PO4
       gamma2(I,26)=gammaHPO4

       gamma2(I,29)=gammaH2OA
       gamma2(I,30)=gammaHOA
       gamma2(I,31)=gammaOA
       gamma2(I,19)=gammak
       gamma2(I,36)=gammaca
       gamma2(I,32)=gammamg
       gamma2(I,37)=gammapo4
       gamma2(I,38)=gammapo4
       
         xx=dabs(fcnnew(ma(6)))
c         if (is.le.3 .or.is.eq.9) 
c     &    write(6,'(A,2I5,5E15.6)')'xx', IS,iii, ma(6),gammah,gammah0
        if(dabs(gammah0/gamma2(Is,6)-1).le. 1D-3.and.xx.le.1D-11 
     &      .and. dabs(xx/ma(6)).le.1D-3) goto 22 
       gammah0=gamma2(is,6)
      enddo



 22   continue
      
         xx=fcnnew(ma(6))
         xx2=xx

         dx= dabs(xx/ma(6))
         if (dabs(xx).ge.1D-10 .and. dx.gt.1D-4) then

               do i=1,1
               call calHNewaa(Ta,ma)
         xx=fcnnew(ma(6))
c         PRINT*, 'calh xx', I, ma(6),XX
      ENDDO


            endif
      
      
       return
       end



      	subroutine calHNew_model(Ta,ma)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50,NSMM=100, npsolid=12)
        real*8 fcnnew
	external  fcnnew

        
        real*8 MA(NP),m(NP),xn(NSMM, np)
        common /xn/xn

	common/suls/ T, M, ppartco2

       
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

       common /kout/xx2
       
      common /output/imode_output,idiff,imode_pH,imode_eq

       common /time/ time, xn_area
       common /test/ Ierr
       real*8 gamma2(NSMM*2,NP)
      common /gamma2/gamma2
c      print*,'is in calhnew',is
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     &     ,gammaH2OA       ,gammapo4
        
        
      if (IS.gt.2*NS) IS =1
      if (ns.lt.1) IS=1
      if (Is.lt.1) IS=1

      ii=0
      DO j=1,NP
c         gamma2(Is,j)=1d0
      enddo
      gammah0 = gamma2(IS,6)

      
       DO iiI=1,10
       
c interat H+
         ii=0

         call calHNewaa(Ta,ma)

         call aw_back_model
     &         (ta,Ma,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         I=IS
       gamma2(I,13) =get_gammaco2(ta,Ma)       
       gamma2(I,6)=gammah
       gamma2(I,12)=gammaNH4
       gamma2(I,17)=gammaCl
       gamma2(I,16)=gammaNa
       gamma2(I,18)=gammaNO3
       gamma2(I,27)=gammas1
       gamma2(I,28)=gammaS2

       gamma2(I,8)=gammaNO3
       gamma2(I,34)=gammaNO3
       
       gamma2(I,1)=aw
       xmplus=ma(6)+ma(19)+ma(12)+ma(16)
       gammaOH =get_gammaOH(xMplus)
       gamma2(I,7)=gammaoH

       gamma2(I,15)=gammaHCO3
       gamma2(I,14)=gammaCO3

       gamma2(I,24)=gammaH3PO4
       gamma2(I,25)=gammaH2PO4
       gamma2(I,26)=gammaHPO4
       gamma2(I,37)=gammapo4
       gamma2(I,38)=gammapo4

       gamma2(I,29)=gammaH2OA
       gamma2(I,30)=gammaHOA
       gamma2(I,31)=gammaOA
       gamma2(I,19)=gammak
       gamma2(I,36)=gammaca
       gamma2(I,32)=gammamg
         xx=dabs(fcnnew(ma(6)))
c         if (is.le.3 .or.is.eq.9) 
c     &    write(6,'(A,2I5,5E15.6)')'xx', IS,iii, ma(6),gammah,gammah0
        if(dabs(gammah0/gamma2(Is,6)-1).le. 1D-3.and.xx.le.1D-11 
     &      .and. dabs(xx/ma(6)).le.1D-3) goto 22 
       gammah0=gamma2(is,6)
      enddo



 22   continue
      
         xx=fcnnew(ma(6))
         xx2=xx

         dx= dabs(xx/ma(6))
         if (dabs(xx).ge.1D-10 .and. dx.gt.1D-4) then

               do i=1,1
               call calHNewaa(Ta,ma)
         xx=fcnnew(ma(6))
c         PRINT*, 'calh xx', I, ma(6),XX
      ENDDO


            endif
      
      
       return
       end
            


c     -----------------------------------------------------------------------------

 
	subroutine dzbren(func, erabs,tol,x1,x2,ITMAX)
c	FUNCTION dzSbren(func,x1 ,x2,tol)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	save fa,fb,a,b,c,fc,xm,d,e
	common /test/ierr
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)	
        ierr=0
        a00=a
        b00=b
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.)) then
           ierr=1
        print*, 'a',a,fa
        print*, 'b',b,fb
        a=a00
        b=b00
	fa=func(a)
	fb=func(b)	
        print*, 'a0',a,fa
        print*, 'b0',b,fb
           
        endif

	c=b
	fc=fb	
	do  iter=1,ITMAX
	if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
	endif
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	 !Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !	Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.b
	e=d
	endif
	else !	Bounds decreasing tOo slowly, use bisection.
	d=xm
	e=d
	endif
	a=b	!Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
	enddo 
c	print*,'dzbrent exceeding maximum iterations'
	x2=b
	return
	END
      


	subroutine dzbrens(func, erabs,tol,x1,x2,ITMAX)
c	FUNCTION dzbren(func,x1 ,x2,tol)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	save fa,fb,a,b,c,fc,xm,d,e
	
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)

c        print*, 'ass',a,fa
c        print*, 'bss',b,fb
        
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.)) then
           
        print*,  'root muss be bracketed for dzbrents'
        stop
          endif
           c=b
	fc=fb	
	do  iter=1,ITMAX
	if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
	endif
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	 !Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !	Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.b
	e=d
	endif
	else !	Bounds decreasing tOo slowly, use bisection.
	d=xm
	e=d
	endif
	a=b	!Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
	enddo 
       
	print*,  'dzbrents exceeding maximum iterations'
        stop
        
	x2=b
	return
	END
      
      



	subroutine dzbren3(func, erabs,tol,x1,x2,ITMAX)
c	FUNCTION dzbren(func,x1 ,x2,tol)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	save fa,fb,a,b,c,fc,xm,d,e
	
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)

        
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.))
     *       then
       print*, 'a3',a,fa
        print*, 'b3',b,fb

         print*,  'root muss be bbren3racketed for zbrent3'
c         pause
         
      endif
      
         c=b
         fc=fb

         
	do  iter=1,ITMAX


       if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
      endif
      
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif

  
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
        
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	 !Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !	Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.b
	e=d
	endif
	else !	Bounds decreasing tOo slowly, use bisection.
	d=xm
	e=d
	endif
	a=b !Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
      enddo
      
	print*, 'dzbrent3 exceeding maximum iterations'
	x2=b
	return
	END
      
      
      


c     --------------------------------------------
c
c     required by calHnew , calculated the H+ and OH- concentration to ensure to charge balnace
c     
c     ------------------------------------------



c     c   Syed Tasleem Hussain, Gul Abad Khan, Muhammad Shabeer 6 points
c      
c     + Chapin 1931 3 data points
      function apH2C2O4_2H2O(T)
        IMPLICIT REAL*8 (A-H,m,O-Z)
       	fcal=  -0.6226  -0.15309E5*(1/t-1/298.15d0)
        apH2C2O4_2H2O= dexp(fcal)

        apH2C2O4_2H2O= dexp(fcal)*1D-9*1013.5 ! hPA
        return
      end
      
            function cal_sigma(T,RH)
        IMPLICIT REAL*8 (A-H,m,O-Z)
        
        RH1=RH
        if (RH.ge.0.99999) RH1=0.99999
        if (RH.le.0.01 ) RH1=0.01
        
        xx= (dlog(1-RH1)*10)-dlog(1-0.9d0)*10

        
        
        
        cal_sigma=72-40*(datan(xx)/3.14159*2+1)/2
        
        
        
              return
              end


      function apna2hpo4_7h2o(T)
      implicit real*8 (a-z)
      
c     aplog= -6.333 + (-5.57515+6.333)/5d0*(t-293.15d0)
c      aplog= -6.333 + (-5.57515+6.333)/5d0*(t-293.15d0)
      aplog= 46.1433-1.53887E4/T
c     aplog= 47.807-1.5551E4/T
      aplog=  0.41777E+02       -0.13673E+05 /t
      apna2hpo4_7h2o=dexp(aplog)
      return
      end
      

c     Dos Santos 2019
      function get_gammaco2(t,M)
      implicit real*8 (a-h,o-z)
c
      common/lam/xlamNa,xlamcl
      
      real*8 m(*)
      
      real*8 ANa(9)
      real*8 ACl(9)
      DO kk=1,9
      ANa(kk)=0
      ACl(kk)=0
      enddo
      Ana(1) =1.694
      Ana(2) = -8.616E-3
      Ana(3) = 1.824E-5
      Ana(4) = -1.361E-8
      Ana(5) = -87.091
cc      acl(1)= .916
c      acl(2)= -5.719E-3
c      acl(3)= 1.351E-5
c      acl(4)= -1.144E-8
      
      xm16=m(16)
      xm17=m(17)
      if (xm16.ge.6d0) xm16=6d0
      if (xm17.ge.6d0) xm17=6d0

      
      etaNaCl=-0.002
      y=0
      y=xm16*xm17*etaNaCl
c     Na+

      xlamNA =ANA(1)+ana(2)*t+ ana(3)*t**2+ana(4)*T**3+ana(5)/T
      xlamCl = 0d0! Acl(1)+acl(2)*t+ acl(3)*t**2+acl(4)*T**3
      y=y+ 2*xm16 * xlamNa + 2*xm17*xlamCl
      y=dexp(y)
      get_gammaco2=y
c      print*, xlamNa
c      print*, xlamCl
c      print*, mc(3), ma(2)
      
      return
      end
      
      
c      call Dh2O_NaCl(T,xm,aw,dl)

       subroutine DH2O_Nacl_suc(T,ML,aw,dl)
      implicit real*8 (a-h,o-z)
      parameter (np=50)
      real*8 M(NP),ml(NP)
      
      real*8 x3(3),y3(3),x1(1),y1(1)
      
      data a0/ 0.21319213/
      data a1/ .13651589D-2/
      data a2/ -0.12191756D-5/
      data b0/ 0.069161945/
      data b1/ -0.27292263D-3/
      data b2/0.20852448D-6/
      data c0/-0.25988855D-2/
      data c1/0.77989227D-5/

      data x3 /0, .4, .74/
      logical ex
      
c       data y3 /-5., -2.5, 0d0/
       data y3 /-3., -2., 0d0/

       data key /0/
       save x3,y3,x11,x22
       
       if (key.eq.0) then
      x11=.8
      x22=1d0
          key=1
      INQUIRE (FILE='b_diff.dat', exist=ex)
      if(ex) then
          open(99,file='b_diff.dat')
          DO I=1,3
             read(99, *, end=98) x3(I),y3(I)
          enddo
          read(99,*) x11
          read(99,*) x22
          
 98       continue
       endif
       
       endif
       



      xm = aw2m(aw)
      
      
      if (xm.ge.6.3d0) xm=6.3d0
 
      a= a0+a1*t+a2*t*t
      b= b0+b1*t+b2*t*t
      c= c0+c1*t

c     do not use c
      c=0d0
      ratio = dexp(a *xm+ b*xm**2)



       aw0=1d0

       dl0= dh2o_func(t)

       
       dl=dl0/ratio

c      m(1)=1000/18d0
c           m(16)= 6d0
c           call aw_back
c     & (t,M,aw1,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      aw1=.76
      
      if (aw.le.aw1) then
         N3=3
         N1=1
         x1(1)=aw
         call intpl(x3,y3,n3,x1,y1,n1)
         dl=dl*dexp(y1(1))
         endif
      
       DL1=dl
c     calculate binury sucrose
c      print*,'aw ', aw, dl1
      
      call   cal_dlaw_suc(T,aw,dl2)      


c    12   0   c
c    13   1   


      ratio = ML(16) /( ML(2)+ml(16))
c      ratio = ML(16)*58.44 /( ML(2)*342+ml(16)*58.44)

      alpha=1d0
      

      
      

      if (ratio.lt.0.99999d0) alpha =dexp ( x11* (1-ratio)**x22)
      xx=alpha* ratio
      if (xx.ge.1) xx=1
      dl = dl2**(1-xx) * dl1**(xx)
c      print*, 'aw,ra', aw, ratio, dl
      
      return
      end
      

          function aw2m(aw)
                   implicit real*8 (a-h,o-z)
                   
                   if (aw.le..74) then
                      aw2m =6.25
                      return
                   endif
         
                   xx=dlog(aw)
                   aw2m =   -0.38080E+02*xx-0.10952E+03*xx*xx
     &                  -0.17598E+03*xx**3
                   

                   return
                   end
      
c     diffusion coefficients in cm-2 of Nacl, NaHCO3, Na2CO3 mixtures
      

      subroutine DH2O_NaCO3_cl(T,aw,ml,dl,dlion)
                   implicit real*8 (a-h,o-z)
                   parameter (np=50,npsolid=12) ! number of species
                   real*8 ml(*),m(NP),b(5)
                   data key /0/
                   data b /0d0,-0.64900E+01,-0.856, +0.55890E+01,0d0/
      logical ex
      save fion,scale, ex,b
      
      
                   if (key.eq.0) then
                   fion=0.
       scale=1
                   key=1
        INQUIRE (FILE='b.fion', exist=ex)
       if (ex)then
       open(99,file='b.fion')
       DO I=1,5
       read(99,*,end=33) ii, b(I)
      enddo
      
       read(99,*,end=33) scale
 
 33    close(99)
       endif

       
      endif
                 
c                   xmcl=aw2m(aw)
c                    m=0
                    DO kk=1,np
                       m(kk)=0d0
                    enddo
                    
                    m(1)=mL(1)
                    m(16)= ml(16)
                    m(17)= ml(16)
                    
c     get Dl for H2O in binary NaCl solution

      call DH2O_Nacl_suc(T,M,aw,dlNacl)
      call DH2O_NaCO3(T,ml,dlCO3)



      
      ratio = ML(17) /( ML(17)+ml(15)+2*ml(14)) ! NaCl mole fraction
c      ratio = ML(16)*58.44 /( ML(2)*342+ml(16)*58.44)
      alpha=1d0
      

         if (ratio.lt.0.99999d0) alpha =dexp ( 0.8* (1-ratio))
         xx=alpha* ratio
        if (xx.ge.1) xx=1

        dl = dlCO3**(1-xx) * dlNacl**(xx)

c     print*, 'aw,ra', aw, ratio, dl
       call cal_dlaw(T,aw,dlI)  !diffusin coefficient of ions
       
       call cal_dlaw_walker_h2o(T,aw,dl0)

c       if (ex) then
c       freduc=(DLI/dl0)**fion   ! less than SLF
c      else
c    1    0.00000E+00    0.00000E+00    0.00000E+00
c    2   -0.64900E+01   -0.64900E+01   -0.64900E+01
c    3   -0.85600E+00   -0.85600E+00   -0.85600E+00
c    4    0.55890E+01    0.10000E+02   -0.10000E+02
         aw1=aw
         if(aw1.le.0.5) aw1=0.5
         xx= 1-aw1
         freduc=  b(2)*xx+b(3)*xx**2 +b(4)*xx**3+b(5)*xx**4
         freduc=dexp(freduc)


      
       dlion= dl * freduc ! less than SLF
       
       dl=dl*scale
       dlion=scale*dlion
                    return
                    end
      
c     viscosity
      
      subroutine DH2O_NaCO3(T,ml,dl)
                   implicit real*8 (a-h,o-z)

                   real*8 ml(*)
                   
                   
                   xmNa= ml(16)
                   aw0=1d0
c                    call  cal_dlaw_citric(t,aw0,dl0)
                   dl0=dh2o_func(T)
                    
       ratio2 = dlog(120d0)/20d0 * xmNA
       ratio1 = dlog(7d0)/10d0 * xmNA
       
       ratio = ratio2 * ml(14)/(ml(14)+ml(15)+1D-10) +
     & ratio1* ml(15)/(ml(14)+ml(15)+1D-10)
       
       DL = dl0/dexp(ratio)
       
       return
       end


      	subroutine calHNewaA(Ta,ma)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=50,NSMM=100)
        real*8 fcnnew
	external  fcnnew

        
        real*8 MA(NP),m(NP),xn(NSMM, np)
        common /xn/xn

	common/suls/ T, M, ppartco2

       
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq

       common /kout/xx2
       

      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3

       common /time/ time, xn_area
       common /test/ Ierr

       Ierr=0
       isupdate=0
       
       iseqco2=0

	t=ta
        DO kk=1,np
           M(kk)= MA(kk)
           enddo

        
        

c     H+ concentrations
        if (imode_ph.eq.2 ) then
        mA(6)=1D-7
        mA(7)=1D-7
c
        x6=1D-7
	xx6= fcnnew(x6)
           do kk=1,np
              ma(kk)=m(kk)
c              print*, kk, ma(kk)
           enddo
c           call aw_back
c     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

           if (imode_NH4NO3.eq.1) then
              ma(12)=ma(4)
              ma(23)=0d0
           endif
           
              
        return
        endif

        
c        print*,'  is',is ,m
        
        if (is.le.0 .or. ma(6).le.1d-19) then
        xmin=1.E-19
 	xmax=1D-7        ! 100% dissociation of acetic acid
        endif

        if (is.ge.1 .and. mA(6).GE.1d-19) then
           x6 = M(6)
           
        
           xmin=x6/2d0
           xmax=x6*2d0

	xx6= fcnnew(x6)

        dd=.1D0


        if (xx6.ge.0d0) xmin=x6
	xx1= fcnnew(xmin)

c     decrease xmin until xx1 >0

        
        if (xx1.le.dd) then
        DO jj=1,100
           xmin=xmin/2
           xx1= fcnnew(xmin)
           if (xx1.ge.dd .or. xmin.le.1D-19) goto 333
          enddo
 333      continue
          endif


	xx2= fcnnew(xmax)
c     increase xmax until xx2 <0
        if (xx2.gt.-1*dd) then
           DO JJ=1,100
           xmax=xmax*2
	xx2= fcnnew(xmax)
        if (xx2.le.-1*dd .or. xmax.ge.40) goto 223
        enddo
 223    continue
        endif


	endif



        
	erabs=0.d0
	errel=0d0
	ITmax=500

        if (is.ge.1 ) then
c           print*, 'xmin 1', xmin,xx1,is
c           print*, 'xmax 1', xmax,xx2,time
c           print*,  'xn6 ', xn(IS,6)

           endif


        if (is.ge.1 .and. xn(Is,6).gt.0d0) then

c           print*, 'xmin', xmin,xx1,is
c           print*, 'xmax', xmax,xx2,time

           else
	xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)
           endif
           xmax1=xmax

           xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)



	if(xx1.le.0d0) then
	xmax=xmin
	goto 2
      endif
      if(xx2.ge.0.)goto 2

c        print*,' calhneq axx1,xx2', xx1,xx2

        call dzbren(fcnnew, erabs,errel,xmin,xmax,ITMAX)
 2      continue
c        print*, xmax
        if (ierr.eq.0)then
           isupdate=1
        else
           isupdate=0
           endif
        M(6)=xmax
        xx2= fcnnew(xmax)
        if (ierr.ne.0) then
           print*,'ierr =' , ierr
            write(39,'(A,i5)') 'ierr =' , ierr
           write(39,'(A,2E14.4,I5)') 'calhnew', xmax,xx2,is

           write(6,'(A,2E14.4,I5)') 'calhnew', xmax,xx2,is
           print*,'xmin,xmax', xmin, xmax1
c           write(39,'(A,2E14.4,I5)') 'xmin,xmax1', xmin, xmax1
           do kk=1,np
              m(kk)=ma(kk)
              enddo
c     M= MA
        xx1= fcnnew(xmin)
        xx2= fcnnew(xmax1)
        print*,'xx1,xx2 ', xx1,xx2
           write(39,'(A,2E14.4,I5)') 'xx1,xx2', xx1,xx2

           DO kk=1,NP
              m(kk)=ma(kk)
           enddo

c        m=MA

        
           DO I=1, NP
c              print*, i, ma(I)
              enddo

c              print*,'is', is

           endif


           DO kk=1,NP
              ma(kk)=m(kk)
           enddo
c     MA=M
c        is=0
        
 999            return
	end  

c     -------------
      function apKCL(T)
      implicit real*8 (a-z)
      apKCL=    0.93705E+02  -0.52159E+04 /T   -0.13019E+02 *dlog(t)
      apKCL=dexp(apKcl)
      
      return
      end

      function apNaH2PO4(T)
      implicit real*8 (a-z)
      dd=      -0.54237E+00  -0.56171E+03/T   
      dd= 0.93159E+01-0.32027E+04 /t
      
      apNaH2PO4=dexp(dd)
      
      return
      end

     
      
      
      function Dh2O_func(T)
      implicit real*8 (a-z)
      dH2O_func = 0.39300E+03 *dexp(-1073.4/(T-90.9d0))
      
c    1    0.39300E+03    0.39300E+03    0.39300E+03
c    2    0.10734E+04    0.97321E+02    0.10742E+05
c    3    0.90900E+02    0.90900E+02    0.90900E+02
      dH2O_func=dH2O_func*1D-5  ! cm2/s
      
      
      return
      end
      
      
c1): K. KRYNIC,Pressure and temperature dependence of self-diffusion in water , 
cFaraday Discussions of the Chemical Society , 66, 1978, 199-207
      
      subroutine matrix_inv(a,b,xb,x,n)
    ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
      implicit real*8 (a-h,o-z)

      integer  i,j,k,l,m,n,irow
       real*8 big,a(100,100),b(100,100),x(100),dum,xb(100)

      !build the identity matrix
       do i = 1,n
          do j = 1,n
             b(i,j) = 0.0
          end do
          b(i,i) = 1.0
       end do
      do i = 1,n ! this is the big loop over all the columns of a(n,n)
      ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
      ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
      big = a(i,i)
      do j = i,n
         if (a(j,i).gt.big) then
            big = a(j,i)
            irow = j
         end if
      end do
! interchange lines i with irow for both a() and b() matrices
      if (big.gt.a(i,i)) then
         do k = 1,n
            dum = a(i,k)        ! matrix a()
            a(i,k) = a(irow,k)
            a(irow,k) = dum
            dum = b(i,k)        ! matrix b()
            b(i,k) = b(irow,k)
            b(irow,k) = dum
         end do
                end if
 ! divide all entries in line i from a(i,j) by the value a(i,i);
 ! same operation for the identity matrix
                dum = a(i,i)
                do j = 1,n
                   a(i,j) = a(i,j)/dum
                   b(i,j) = b(i,j)/dum
                end do
 ! make zero all entries in the column a(j,i); same operation for indent()
                do j = i+1,n
                   dum = a(j,i)
                   do k = 1,n
                      a(j,k) = a(j,k) - dum*a(i,k)
                      b(j,k) = b(j,k) - dum*b(i,k)
                   end do
                end do
             end do

 ! substract appropiate multiple of row j from row j-1
             do i = 1,n-1
       do j = i+1,n
       dum = a(i,j)
      do l = 1,n
      a(i,l) = a(i,l)-dum*a(j,l)
      b(i,l) = b(i,l)-dum*b(j,l)
      end do
      end do
      end do

       DO I=1,n
      xx=0
      DO J=1,n
         xx=xx+ b(I,J)* xb(j)

      enddo
      x(I)= xx
      
c      print*,i,xx,b(I)
      enddo
      

         end


      function ap_solid(Ip, T)
	implicit real*8 (a-h,o-z)
        ap_solid=0d0
        
        if (ip.eq.1) then
           ap_solid=apNH4HC2O4_05H2O(T)
           return
        endif

        if (ip.eq.2) then
           ap_solid=apNH42C2O4(T)
           return
        endif
        if (ip.eq.3) then
           ap_solid=apNaHC2O4(T)
           return
        endif
        if (ip.eq.4) then
           ap_solid=apNa2C2O4(T)
           return
        endif
        if (ip.eq.5) then
           ap_solid=apH2C2O4_2H2O(T)
           return
        endif
        
        if (ip.eq.6) then
           ap_solid=apNa2HPO4_7H2O(T)
           return
        endif

                if (ip.eq.7) then
           ap_solid=apKCl(T)
           return
        endif
        if (ip.eq.8) then
           ap_solid=apNaH2PO4(T)
           return
        endif

        if (ip.eq.9) then
           ap_solid=apNaCl(T)
           return
        endif

        if (ip.eq.10) then
           ap_solid=apCaCo3(T)
           return
        endif

        if (ip.eq.11) then
           ap_solid=apCaHPO4(T)
           return
        endif
                if (ip.eq.12) then
           ap_solid=apmgCo3(T)
           return
        endif
        
        
        return
        end

      function ss_solid(Ip, T,ML)

      implicit real*8 (a-h,o-z)
      parameter (npsolid=12,np=50,NSMM=100)
      integer NSP(npsolid)
      integer NSP_index(npsolid,10)
      real*8 xNSP_nv(npsolid,10)
      COMMON /nsp/ NSP, NSP_INDEX
      COMMON /xnsp/ xNSP_NV
      real*8 gamma2(NSMM*2,NP)
      common /gamma2/gamma2
      real*8 ml(np)
      
       common /isco2/iseqco2,is,ns,isupdate,iseqNH3,iseq
           ap=1d0

      if(IP.ge.3) then
           call calHNew(T,mL)
           call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           DO I=1, nsp(IP)
              J=nsp_index(IP,i)
c              print*, j, xnsp_nv(IP,i),is
              
              if (j.eq.1) then
                 ap=ap* aw ** xnsp_nv(IP,i)
              else
                 ap= ap*(ml(j)*gamma2(is,j))** xnsp_nv(IP,i)
              endif
c              print*,  ml(gotoj),gamma2(is,j)
c              if (IP.eq.12) write(6,'(A,2I5,5E15.6)')
c     & 'AP12', I, J, ml(j),ap              
           enddo
              
        endif
      if(IP.lt.3) then
       call vapnew(T,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,poa,pl)
           call aw_back
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           if (IP.eq.1)
     &             ap=pNH3*poA* dsqrt(aw)         
           if (IP.eq.2)
     &             ap=pNH3*poA* pnh3
       
      endif
      
      
      ss_solid=  ap/ ap_solid(ip,T)
        return
        end


      function apCaCO3(T) ! handbook CRC
        IMPLICIT REAL*8 (A-H,m,O-Z)
        apcaco3 = 3.36E-9
        
c     increase by a factor 6 so that the MEM and SLF recipes are in equilibrium 
        apcaco3 =         apcaco3* 6
        

        return
      end
      


      function apCaHPO4(T) ! Chow LC 2001 solubity of calclium phosphate
        IMPLICIT REAL*8 (A-H,m,O-Z)
        apCaHPO4 = 10d0**(-6.9)
        return
      end

C Pascale Bénézeth, Giuseppe D. Saldi, Jean-Louis Dandurand, Jacques Schott,
C Experimental determination of the solubility product of magnesite at 50 to 200°C,
C Chemical Geology, Volume 286, Issues 12,

      function apmgCO3(T)       
        IMPLICIT REAL*8 (A-H,m,O-Z)
        apmgco3 = 1.585E-8

        return
      end
      

      

 
c      
      subroutine cal_flux(time,x,f,dtime,NS)
       implicit real*8 (a-h,m, o-z)
       parameter (np=50,nsmm=100,npsolid=12)
        real*8 mm(NP),mv(NP),ml(NP),ml0(NP)

        common /buffer/ ph3eq,istrace

       
       real*8 vshell(nsmm),phshell(NSMM),
     & awshell(nsmm),xhshell(NSMM),ml6shell(NSMM),sshell(NSMM,npsolid)
    
       common /diagnose/ Idiagnose
       
      real*8 dl_factor2(2,np),dl_factor(NP)
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin
      common /timeeq/ timeeq,timeeq0

        common /sigma/ sigma 


       common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &      partvaphcl, partvapOA,partvaplac
                    
       real*8 x(*),xn(NSMM,np),f(NSMM,np+npsolid), xm2(nsmm,np)
       
       real*8 w1(NSMM),vap(np)
       
      real*8 g(NSMM,np),c(NSMM,np),a2(NSMM,np),cm2(NSMM,np)

      real*8 velarr(nsmm),xmbufferI(nsmm),aw2(NSMM)
      
      common /cm2/velarr,cm2
      common /pi/ pi
c      in izc(nsmm,np)
      integer izc(np)

      common /M/ MM,mv,izc

      common /awshell/ awshell,xhshell,ml6shell,
     & sshell,phshell,vshell



        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     & , gammah2OA

       real*8 gamma2(NSMM*2,NP),xmtracemax(NP)
       common /gamma2/gamma2


      common /solid/xnsolid
      common /time/timmm, xn_area


      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
      
      
      common /enhance/radius,r1,enh_factor

      common /Ienhance/Ienh,iscenter,isAHO

      common /xn/xn
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqAA

      common /iskinetic/iskinetic, istrue(NP)

      
       common /venti/ venti,fheat_11
       common /vap/vap
      integer NSP(npsolid)
      integer NSP_index(npsolid,10)
      real*8 xNSP_nv(npsolid,10)
      COMMON /nsp/ NSP, NSP_INDEX
      COMMON /xnsp/ xNSP_NV
      integer IDIS(Np)
      
      data keys/0/
      save IDIS,key0
       common /dxmin/ dxmin
       common /xvion/ xvion,dl_ion,dl_ions(NSMM)
       
       common /ishell/ Ishelleq,imode_shell

c     default:        NS+1
c     time > timeeq:  1
c     ismode_NH4NO3: calculated based on thickness and diffusion coefficient  of H2O
       real*8         dlshellH2o(NSMM)
       
       common /fvap/ alpha1,alpha0,ppvap,fvap_factor,vapratio
       
       
      parameter (Ntimes=80000)
      real*8 timetr0(ntimes),rhtr0(ntimes),gash2otr0(ntimes)
      common /ntr0/ Ntr0
      common /ntr/  timetr0,rhtr0

      dl_ion=1D1
      
       Ishelleq=NS+1
       
       dh2oliq=1000
       
      if (key0.eq.0) then
         key0=1
         DO J=1,NP
            iDIS(j)=0
         enddo
c     set disociate species to 1
         idis(6)=1 !H2O
         idis(7)=1

         idis(8)=1 ! acetic acid
         idis(9)=1

         idis(12)=1 !NH3
         idis(23)=1
c
         idis(13)=1
         idis(14)=1
         idis(15)=1

         if( iskinetic.eq.1) then
         idis(13)=0
         idis(14)=0
         idis(15)=0

         endif
         
         
         idis(24)=1
         idis(25)=1
         idis(26)=1

         idis(27)=1
         idis(28)=1

         idis(29)=1
         idis(30)=1
         idis(31)=1

         idis(33)=1
         idis(34)=1
         
      endif


      nss=ns
       npliq=Np-npsolid
       n1=1
       
c     relocat solids the maximal superdaturation to inner 
       if (iscenter.ne.1) then
       ff=2
       DO KS=1,npsolid
          J=npliq + ks
          DO I=3,ns
             if (xn(i,j) .gt.1D-50) then

                DO II=i-1,1,-1
                   if (xn(ii,j).ge.1D-50) goto 459
                   
                   if ( sshell(ii,ks).ge.ff*sshell(I,ks) .and.
     &                   sshell(Ii,ks).ge.1.5d0)then
                      iii= II
                      if(ii.eq.1) iii=2
                      print*,'relocat solid ks, from i, i1',
     &                     KS, I,iii
                      xn(iii,j)=xn(iii,j)+ xn(I,j)
                       xn(I,j)=0d0
                       goto 459
                    endif

                 enddo
                
             endif
          enddo
          
 459   continue
       enddo
       
c     relocate to outer shells with higher s
       DO KS=1,npsolid
          J=npliq + ks
          DO I=1,ns-1
             if (xn(i,j) .gt.1D-50) then

                DO II=i+1,NS-1
                   if (xn(ii,j).ge.1D-50) goto 469
                   ssI=max(sshell(I,ks),sshell(I-1,ks))
                   if (sshell(ii,ks).ge.ff*ssI .and.
     &                   sshell(Ii,ks).ge.1.5d0)then
                      print*,'relocat solid ks, from i, i1',
     &                     KS, I,ii
                      xn(ii,j)=xn(ii,j)+ xn(I,j)
                       xn(I,j)=0d0
                       goto 469
                    endif
                enddo
                
             endif
          enddo
          
 469   continue
       enddo
       


      endif
      


       
       
       

      
c     calculate the thickness and concentration in each shell


       vols=0d0 ! cm3
       mL(1)=1000/mm(1)
       DO I=1,Ns

          do j=1,NP
c             xm2(i,j)= xn(I,j)/xn(I,1)*mL(1)
          enddo

           if (time.ge.timeeq) ta=t
              ff= 1000d0 /(MM(1)*xn(I,1))
               ML(1)=1000d0/MM(1)
                 DO J=2,NP
                  ML(J)= ff* xn(I,j)
                  enddo
           IS=I
        call calHNew(Ta,mL)
        call aw_back(ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
        call cal_MV(aw,ta,mv)
                    
          vol= mv(1)*xn(I,1)+mv(2)*xn(I,2)

             DO J =6,20
             vol=vol+xn(I,j)*mv(j)
             enddo
c     exclude Titers

             DO j=23,np
               vol=vol+xn(I,j)*mv(j)
               enddo

              vols=vols+vol



      DO jj=1,np
      Ml0(jj)=ml(jj)
      enddo
        
c        Ml0=ML
c     check aw
c         sigma = 72
         dx= x(NS+1)-x(NS)

                
         gamma2(I,19)=gammaK
         gamma2(I,32)=gammaMg
         gamma2(I,36)=gammaCa
         
         gamma2(I,6)=gammah
       gamma2(I,12)=gammaNH4
c       gammaHCO3=gamma2(I,15)
c       gammaCO3=gamma2(I,14)
       gamma2(I,17)=gammaCl
       gamma2(I,16)=gammaNa
       gamma2(I,18)=gammaNO3
       gamma2(I,27)=gammas1
       gamma2(I,28)=gammaS2
c       print*, 'g',       gamma2(I,26),gammaHPO4

       gamma2(I,25)=gammaH2PO4
       gamma2(I,26)=gammaHPO4
       gamma2(I,29)=gammaH2OA

       gamma2(I,30)=gammaHOA
       gamma2(I,31)=gammaOA

       gamma2(I,1)=aw
       xmplus=(xn(I,6)+xn(I,12)+xn(I,16))/xn(I,1)*55.51
       gammaOH =get_gammaOH(xMplus)
       gamma2(I,7)=gammaoH
       

c        gammaOH =get_gammaOH(xMplus)
c        xmplus= ml(6)+ml(16)+ ml(12)+ml(19)
         
        a2(I,1)= aw
        DO j=2,np-npsolid
        a2(i,j)= ML(J)*gamma2(i,j)
        enddo
        


         awshell(I)=aw
         Amisch1= ML(16)*ml(17)*gammaCl*gammaNa


         xhshell(I)=gammah
         ml6shell(I)=ML(6)
         phshell(I)= dlog(ML(6)*gammah)/dlog(.1d0)
c         print*,'ph', i,phshell(I)
         
        call vapnew(Ta,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,pOA,pl)

        DO IP=1,npsolid
           sshell(I,ip) = ss_solid(ip,Ta,ml)
        enddo
        ss= apna2hpo4_7h2o(T)

c        write(6,'(A,I5,10E15.6)') 'ss6',I,ss
c     & aw**7* ml(16)**2 *gamma2(I,16)**2 * ml(26) *gamma2(i,26),
c     & aw, ml(16) ,gamma2(I,16) , ml(26) ,gamma2(i,26)

          if (I.eq.NS) then
         vap(1)=aw*vwater(T) ! take the air T head considered via Pruppacher
         vap(3) = pacetic
         vap(4) = pNH3
         vap(5) = pCO2
             vap(18) = pHNO3
          vap(17) = pHCL
          vap(29) = pOA
          vap(33) = pl
          endif
       
         ff= xn(I,1)/ML(1)

c     recalculates equilibrium species
         DO j=6,15
            xn(I,j) = ff*ML(J)
         enddo

         DO J=23,np
            xn(I,j) = ff*ML(j)
         enddo
c     liquid volume

         vshell(I)=0

         DO kk=1,np-npsolid
         if (kk.ge. 23 .or. kk.le. 20) 
     & vshell(I)= vshell(I) + mv(kk)*xn(I,kk)
         enddo


         vcore= 4*pi/3* x(1)**3
         xx = 3d0/4d0/pi* (vols+vcore)
         
              x(I+1) = (xx)**(1/3d0)
c              print*, i,x(I+1)
              DO J=1,NP
              c(I,j)= xn(I,j)/vshell(I) ! mol / cm3 
              enddo              
           enddo
c           stop
           
c     important need for the enhancement factor 
           radius=x(NS+1)
           vv=0
           DO kk=np-npsolid+1,np
              vv=vv+ mv(kk)*xn(1,kk)
              enddo
           r1= (vv/4d0/pi*3d0)**(1/3d0)

           DTIMEMIN= 10
           
          xnsolid=0d0

          DO I=1,NS
             do J0=np-npsolid+1,np
          xnsolid=          xnsolid+xn(I,j0)
             enddo
             enddo



c     NH42(OA) oxalate  crytal
c     kk=NP: NaCl
              
              DO KK=Np-npsolid+1, NP
                 ISAHO= kk-np+npsolid
c     iss= solid1
          DO I=1,NS
c             print*, I, xn(I,kk)
                 f(I,kk)   = 0d0  !flux to solid from inner shell
                 f(I,np+isaho) =0d0  !flux from outer shell to solid
             
              if (xn(I,kk).gt.1D-50 ) then

                 IS=I

c     NH42(OA) oxalate  crytal                 
                 call cal_AO_flux (NS, I,ta, x,ff1,ff2)
c               if (dabs(ff1)+dabs(ff2).gt.0d0 .and. ISAHO.eq.11) then
c                  ff=ml(1)/xn(I,1)
c                  write(6,'(A,2I5,15E15.6)') 'FF2 ',I,isaho,ff2,ff1
c     &  , sshell(i,ISAHO)

                 f(I,kk)   = ff1  !flux to solid from inner shell
                 f(I,np+isaho) = ff2  !flux from outer shell to solid
                 if (idiagnose.ge.4) then
                 print*,'fluxcal', I, ISAHO, ff1,ff2
              endif
                 ml(kk)= ml(1)*xn(I,kk)/xn(I,1)



c     determine gas phase change



              endif

              enddo
           enddo



c             print*, 'befor solid'

c              DO I=2,NS
c                 DO J=1,np
c                    izc(I,j)=izc0(J)
c                 enddo
c              enddo
c     calculate the mean concentration at th boundary

              DO I=2,NS
c     do extra equilibrium calculation when pH step > 0.1

                 DO J=2,NP
c            ml(J)= (xn(I-1,j)/xn(I-1,1)+xn(I,j)/xn(I,1))/2*ml(1)
                 enddo
                 IS=I+NS
c                 call calhnew (ta,ml)
c           call aw_back
c     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

      aw2(I) =(awshell(I)+awshell(I-1))/2
c                 
c                 volm =  mv(2)*ml(2) + mv(1)*ml(1)
c                 DO J=6,np-npsolid
c                    volm=volm+ ml(j)*mv(j)
c                 enddo


                 
                 DO j=1,NP-npsolid
        if (xn(I-1,j).gt.0d0.and.xn(I,j).gt.0d0.and. idis(j).eq.2) then
c      take geometric mean for CO2 and HCO3-
         cm2(i,j)=dsqrt(xn(I-1,j)/vshell(I-1)*(xn(I,j)/vshell(I)))
                     else
          cm2(i,j)=0.5*(xn(I-1,j)/vshell(I-1)+(xn(I,j)/vshell(I)))
                       endif

                       if (imode_ph.eq.2) then

c        if (xn(I-1,j).gt.0d0.and.xn(I,j).gt.0d0) then
c         cm2(i,j)=dsqrt(xn(I-1,j)/vshell(I-1)*(xn(I,j)/vshell(I)))
c      else
c          cm2(i,j)=0.5*(xn(I-1,j)/vshell(I-1)+(xn(I,j)/vshell(I)))
c        endif
                endif

                enddo
             enddo

              
              
              
              
        Ienh=0

        DO J=NP-npsolid,1,-1
c      set the flux at radius = 0 to zero when no solid

           f(1,J)=0
           
           if (dl_factor2(1,J).gt.0 .and. istrue(j).eq.1) then

c     calculates aw
c     Special treatment for solid NaCl for f(1,16) Na+ and f(1,17) Na+ and Cl-
c     f(1,16) = f(1,17) Na+ and Cl- has the same flux
c     D16 * (C1 -c0) / dx
              Ienh=0
                             
       DO I=2,NS

       thickness = (x(I+1) -  x(I-1))/2
c     calculates the gradient in concentration
c
          
          g(i,J)= (c(I,J)-c(I-1,J))/thickness


c          cmm= (c(I,j)+c(I-1,j))/2
          
          cmm= cm2(I,j)


          if (a2(I,j)  .gt.1D-30 .and. a2(I-1,j)  .gt.1D-30) then
c     & (Idiff.le.1 .or. idiff.ge.5) ) then


c     takes the activity/aw as the driving force
c          g(i,J)= cmM/thickness *
c     & (a2(I,J)-a2(I-1,j))/(a2(I,j)+a2(I-1,j))

             g(i,J)= cmM/thickness *
     & dlog( a2(I,J)/a2(I-1,j))

       endif

c     calculates the diffusion coefficients
c     takes the mean value of shell I and I-1
c       aw=(awshell(I)+awshell(I-1))/2 ! take the mean value of the neighboring shells
         dl_factor(J)=dl_factor2(1,J)  
              aw=aw2(I)
c     take the mean value
              DO k=2,np
              ml0(k)=ml(1)/2d0*(xn(I-1,k)/xn(I-1,1)+xn(I,k)/xn(I,1))
              enddo
              DO k=1,NP
c                 print*,k,ml0(k)
              enddo
              
         call caldl(ta,aw,ml0,x,J,dl)
c     for water takes modified r data 
         
         d= dl*dl_factor(J)     ! scaling facotr for species J, 
             
             if (J.eq.1 .and. x(I+1).le. 0.8 *x(NS+1)
     & .and. dh2oliq.gt.d)             dh2oliq=d
c     calculate the flux rate from shell i+1 to i : the direction
c      if (J.eq.1) write(50,'(1F15.4, I5,4E15.6)') time,I, aw,d,
c     & dl_factor(J)       

             
          f(I,j)= g(I,j) *d* 4*pi* x(I)**2 ! the missing minus sign changes the diffusion direction

c     reduction due to solid
          vv= 4*pi/3 *( x(I+1)**3-x(I-1)**3)

          vsolid = mv(NP)*(xn(I,NP) +xn(I-1,NP))+
     &  mv(NP-1)*(xn(I,NP-1) +xn(I-1,NP-1))

          f(i,j) = (vv-vsolid)/vv* f(i,j)


c     for acetic acid take the sum of 8,9 10
                if (dl_factor2(1,8).gt.0 .and. j.eq.3) then
 		f(I,3)=f(I,8)+f(I,9)+f(I,10) !total acetic acid flux
                endif


                
                
c     takes the sum of NH4CH3COO, NH4+, NH3
                if (dl_factor2(1,12).gt.0 .and. j.eq.4) then
 		f(I,4)=f(I,10)+f(I,12) +f(I,23) ! total NH3 flux
                endif

                if (dl_factor2(1,15).gt.0 .and. j.eq.5) then
 		f(I,5)=f(I,13)+f(I,14) +f(I,15) ! total Co2 flux
                endif
                
         if (time.ge.timeeq.and.j.eq.1)f(I,1)=0d0


       enddo              !I=2,Ns


       I=NS
       

c     the gas phase deposition has to be calculated extra.
c     calculates the gas phase deposition! no gas phase depletion is considered, constant gas phase is assumed here
       
      DH2O0=venti*0.211*1013.d0/PRESS*(T/273.15d0)**1.94 ! Pruppacher + Klett

      rad= x(Ns+1)
      alpha=1
      if (j.eq.4 .or.j.eq.18 .and. imode_NH4NO3.eq.1) then
         xmaa=(mm(12)*xn(I,12)+mm(18)*xn(I,18))
         
         ratioan= xmaa/(xn(I,2)*mm(2)+xmaa )
c         alpha= alpha0 + alpha1 * ratiom
        alpha= alpha0+ (alpha1-alpha0)* ratioAN
         
      endif
      
       velocity=(8.d0*8.314D7*T/PI/MM(j))**0.5
        Dh2o =DH2O0/(1d0 +
     +        (4.d0*DH2O0/RAD/alpha/velocity))
      vair = 18.5*.8+ 16.3*.2
      xmd = 2/(1/mm(j)+ 1/28.97)

c     Equation 10 of Langenberg 2020
        DGas1 =DH2O0* dsqrt(MM(1)/MM(J))
        pp=press/1013.5
        
      
      Dgas=venti*0.00143*
     *   T**1.75/dsqrt(xmd)/(vair**(1/3d0) + mv(j)**(1/3d0))**2/pp
c      print*, ' j Dgas', j, dgas, dgas1
      
        dgas =Dgas/(1d0 +
     +        (4.d0*Dgas/RAD/alpha/velocity))
 
         xkelvin = dexp( 2* sigma * 18d0 /(8.314E7*T*rad) )

         f(NS+1,j)=0d0

         
         if (J.eq.1) then
         xkelvinH2O = dexp( 2* sigma * 18d0 /(8.314E7*T*rad) )
         xngaspartial = rh*vwater(T)*1D-4 /8.314/T ! mol/cm3air
         ph2ovap = vap(1) *xkelvin
         xngasvap = ph2ovap*1D-4 /8.314/Ta ! mol/cm3air
         delta = ( xngaspartial-xngasvap) 
cccccccccccccccccccccc

      vapheat=597.3*(273.15/t)**(0.167+3.67d-4*t) !cal/g
      thermcond0=((5.69+0.017*(t-273.15))*1.d-5) ! cal/s/cm/K
      velocity=(8.d0*8.314D7*T/PI/28.97d0)**0.5
       cair =1 /4.183 ! cal/K/g

       rhoair = press/1000  * 1E-1 /8.314d0/T*28.97d0 !  g/cm3

c 
      thermcond=thermcond0/(1.d0+4.d0*thermcond0/
     +                     (RAD*velocity*Rhoair*Cair))
      
      fheat=fheat_11*18d0*vapheat/ thermcond /T *
     &(vapheat*18d0/1.987d0/t-1)*
     &  vwater(t) *100/8.314/T*1E-6 ! mol/cm3
     & *dgas 

      f(NS+1,j) = 4*pi*rad* DGAS  *delta/(1+ fheat)
c      if    (dx.le.4.0d0*dxmin.and. time.ge. .1*(x(NS+1)*1D4)**2)
c     &       f(NS+1,1)=f(NS,1)

      if (time.ge.timeeq) f(NS+1,j)=0d0
       dx=x(NS+1)-x(NS)
       
         
  
      Ta=t + fheat_11*vapheat*mm(1) * f(NS+1,1) /4/pi/rad/thermcond
ccccccccccccccccccccc

       if (imode_nH4NO3.ge.1) ta=t
      
      xx= x(NS+1)-x(NS)
c     set last shell to equilibrium for H2O

         
c     print*, dgas, ph2o0, rad, xngasvap,t,xkelvin
      endif

c     acetic acid deposition
         IF (J.eq.3                     ) then
            
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvap3*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(3)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air
         dgastang= 94d0/760 *1013d0/press *(T/296)**1.75
         dgas=dgastang
         dgas =DH2O0/(1d0 +
     +        (4.d0*Dgas/RAD/alpha/velocity))

         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)
         dx=x(NS+1)-x(NS)
c         if (dx.le.4.0*dxmin.and.time.ge. .1*(x(NS+1)*1D4)**2 )
c     & f(NS+1,1)=f(NS,1)
c              if (time.ge.timeeq)f(NS+1,j)=0d0

         
      endif
      if (imode_NH4NO3.eq.1 .and. (j.eq.4 .or. j.eq.18)) then
         pp= ppvap* 1D-18 *press**2 !hpa**2
         
         partvapHNO3=  (pp/vapratio)**0.5d0
         partvap4=  (pp/vapratio)**0.5d0 *vapratio
         
         
         pp=vap(4)*vap(18)
         
         if (pp.gt.0d0) then
         dgastang= 87d0/760 *1013d0/press *(T/296)**1.75
       velocity=(8.d0*8.314D7*T/PI/MM(18))**0.5

         dgastang =dgastang/(1d0 +
     +        (4.d0*dgastang/RAD/alpha/velocity))

         dgastangNH3= 176d0/760 *1013d0/press *(T/296)**1.75

         velocity=(8.d0*8.314D7*T/PI/MM(4))**0.5
         dgastangNH3 =dgastangNH3/(1d0 +
     +        (4.d0*dgastangNH3/RAD/alpha/velocity))
         fD= dgastang/dgastangNH3
         bb= partvap4 - fD * partvapHNO3
         vap(18) = 1/2d0/fD * ( -bb + dsqrt( bb*bb + 4*pp* fD))
         
         vap(4) = pp/vap(18) !dsqrt(pp)  /dsqrt (dgastangNH3/dgastang  )

c         print*,'part',  time, partvap4,partvapHNO3
c         print*,'fD', fD
c         print*, 'delta P',vap(4)-partvap4, fd*(vap(18)-partvapHNO3)
c               print*, ppvap*1D-18 *press**2, pp
c         stop
         
c     ppvap ppb**2
         

c     calculate partial press from ppavp
         

      endif
      endif
      

      
c NH3
      IF (J.eq.4 ) then
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )
         xngaspartial = partvap4*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(4)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air
         dgastang= 176d0/760 *1013d0/press *(T/296)**1.75
         dd=dgastang

         dgastang =dgastang/(1d0 +
     +        (4.d0*dgastang/RAD/alpha/velocity))

c         print*, 'NH3 ',rad, dd/dgastang
c         stop
         
         dgas=dgastang
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)

         
      endif
c      if (j.eq.4) print*,'f', f(NS+1,j)

      IF (J.eq.555 ) then
c     for CO2 steady state
c     not used
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

          dx=x(NS+1)-x(NS)

       xngaspartial = partvapco2*1D-4 /8.314/T ! mol/cm3air
       xngasvap = vap(5)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air
         
         dgastang= 0.16 *1013d0/press *(T/293.15)**1.75 ! engineering tool box
c         print*, j,dgas,dgastang
         dgastang =dgastang/(1d0 +
     +        (4.d0*dgastang/RAD/alpha/velocity))
         dgas=dgastang
        f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)

        
      endif

       IF (J.eq.18 ) then !HNO3
          xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvapHNO3*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(18)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air

         dgastang= 87d0/760 *1013d0/press *(T/296)**1.75
c         print*, j,dgas,dgastang

         dd= dgastang
         dgastang =dgastang/(1d0 +
     +        (4.d0*dgastang/RAD/alpha/velocity))
c         print*, 'NO3 ', rad, dd/dgastang

         
         dgas=dgastang
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)         
c         print*, '18', f(NS+1,j),f(NS+1,4)
       endif

      
      IF (J.eq.17 ) then         !HCL
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvapHcl*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(17)*1D-4 /8.314/T*xkelvin ! mol/cm3air
          dgastang= 118d0/760 *1013d0/press *(T/296)**1.75
         dgastang =dgastang/(1d0 +
     +        (4.d0*dgastang/RAD/alpha/velocity))

          dgas=dgastang
          
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)
         if (xngaspartial .lt. -1D-40) f(NS+1,j)=0

      endif

      IF (J.eq.29 ) then         !Oxalic acid
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvapOA*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(29)*1D-4 /8.314/T*xkelvin ! mol/cm3air


         
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)
         if (xngaspartial .lt. -1D-40) f(NS+1,j)=0



      endif


      IF (J.eq.33 ) then         !lactif
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvaplac*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(33)*1D-4 /8.314/T*xkelvin ! mol/cm3air
          
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)
         if (xngaspartial .lt. -1D-40) f(NS+1,j)=0



      endif

      
      

         endif
         
      ENDDO !J 


c     flux adjustment for charge balance of liquid phase diffusion
c     f+ + c*n+ = f- - c*n-
c     takes only the species with DL_factor > 0 

c     for trace species reduce flux, only if istrace.eq.1 
c      print*,'ff', time,timeeq,istrace
      if (time.ge.timeeq .and. istrace.gt.1) then
         factor=0.1
         DO J=3,np-npsolid

            xmtracemax(j)=0d0
            DO I=1,ns
               xmm= xn(I,j)/xn(I,1)*ml(1)
               if (xmm.gt. xmtracemax(j)) xmtracemax(j)=xmm
            enddo
         enddo
c     acetic acid
c         print*,'xm3', xmtracemax(3)
         
         if (xmtracemax(3).le.1D-4.and.xmtracemax(4).gt.0d0) then
         J=3
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=8
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=9
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=10
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
          endif


         if (xmtracemax(4).le.1D-4 .and.xmtracemax(4).gt.0d0) then
         J=4
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=23
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=10
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=12
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo

      endif

      J=5
      if (iskinetic.eq.0 .and. xmtracemax(5).le.1D-4 .and.
     &      xmtracemax(j).gt.0d0) then

         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=13
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=14
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
         J=15
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo

      endif

      if (iskinetic.eq.1 .and. xmtracemax(14).le.1D-4) then
         J=14
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
      endif

      j=32
      if ( xmtracemax(j).le.1D-4 .and. xmtracemax(j).gt.0d0) then
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
      endif

      j=36
      if ( xmtracemax(j).le.1D-4 .and. xmtracemax(j).gt.0d0) then
         DO I=1,NS+1
           f(I,j) = f(I,j) /istrace
         enddo
      endif

      j=24
      xmm=xmtracemax(j)+xmtracemax(j+1)+xmtracemax(j+2)
      if ( xmm.le.1D-4 .and. xmm.gt.0d0) then
         DO JJ=j,j+2
         DO I=1,NS+1
           f(I,jj) = f(I,jj) /istrace
         enddo
         enddo
      endif
      j=29
      xmm=xmtracemax(j)+xmtracemax(j+1)+xmtracemax(j+2)
      if ( xmm.le.1D-4 .and. xmm.gt.0d0) then
         DO JJ=j,j+2
         DO I=1,NS+1
           f(I,jj) = f(I,jj) /istrace
         enddo
         enddo
      endif

      
      j=27
      xmm=xmtracemax(j)+xmtracemax(j+1)
      if ( xmm.le.1D-4 .and. xmm.gt.0d0) then
         DO JJ=j,j+1
         DO I=1,NS+1
           f(I,jj) = f(I,jj) /istrace
         enddo
         enddo
      endif

      
      j=33
      xmm=xmtracemax(j)+xmtracemax(j+1)
      if ( xmm.le.1D-4 .and.xmm.gt.0d0) then
         DO JJ=j,j+1
         DO I=1,NS+1
           f(I,jj) = f(I,jj) /istrace
         enddo
         enddo
      endif
      
      

      
         
      endif
      
      
      if (Idiff.gt.0 .and. Idiff .le.4 ) goto 234
      if (imode_NH4NO3.eq.1) goto 234  ! no charge balance
      
      velarr(1)=0
      velarr(NS)=0
      
      DO I=2,Ns

         ff=0
         czz=0
 
        DO J=1,NP-NPsolid
        if ((izc(j)).ne.0) then
          ff=ff+ izc(J)* f(i,j)
          czz=czz+cm2(i,j)* izc(J)*izc(J)
          endif
         enddo





C         VEL = V/ DL_FACTOR

         vel= -ff/czz
c     reduces the residue

         ff=0d0
         Vel12 = 0d0! vel
C     c         print*,vel
         velarr(I-1)= vel
         
             DO J=1,np-npsolid
                if ((izc(j)).ne.0) then
               f(i,j)=f(i,j) + izc(J)*vel*cm2(i,j)
            endif

            enddo




            ENDDO


      
 234  continue
c     set Acetic acid to equilibrium
      

c       izc(14)=-1 ! treat 
c       izc(13)=0 ! treat 
         iseqAA=0
         ff=(x(NS+1)*1D4)**2
         if (phshell(NS).le.3) ff=0.01*(x(NS+1)*1D4)**2
       if (phshell(NS).le.2) ff=.001*(x(NS+1)*1D4)**2
         
         if (phshell(NS).le.pH3eq .and. partvap3.gt.0d0
     &        .and.time.ge.ff .and.
     &        xn(Ns,3)*ml(1)/xn(ns,1).le.1D-3) then

c         xka=1.74D-5            ! dissociation constant of acetic acid !
c     Harned and Ehlers 1933
       xx= -1500.65/Ta-6.50923 * dlog(Ta)/dlog(10d0)-0.0076792* Ta
     & +18.67257d0
         xka= 10.d0**(xx)

         Xha =4000 * dexp(6200*(1/ta-1/298.15d0)) ! Henrys law acetic acid M/bar
c     Sander 2015a
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*Ta*x(NS+1)) )

         xm9 = partvap3/1013.5 *xha/xkelvin 
         gammaa=1d0
         xm8 = xka*xm9/(ml6shell(NS)*xhshell(NS))/gamma2(NS,18)
         xm3 = xm8+xm9
         xn(NS,3) =xm3/ML(1)* xn(NS,1)
cc
         f(NS,3) = f(ns,8)+f(ns,9)+f(ns,10)
         f(NS+1,3) =f(NS,3)
         iseqAA=1
c     make equilibrium for HNO3
         endif
      
      iseqNH3=0
c     do NH3 only when alkaline and minor
      
      if (phshell(NS).gt.7d0 .and. partvap4.gt.0d0 
     &     .and.time.ge.timeeq/2.and.
     & xn(NS,4)/xn(ns,1)*ml(1).le.1D-3 .and. imode_NH4no3.ne.1) then

      iseqNH3=1

              xkNH4bb=xknh4b(Ta)*gammaH/gamma2(NS,12)

         xkelvin = dexp( 2* sigma * MV(4) /(8.314E7*Ta*x(NS+1)) )

c         p=m0(12)/m0(6)/ xknh3(T)*1013.5*gammaNH4/gammaH
         xm12= partvap4/1013.15 * ml6shell(NS) *
     &    gamma2(NS,6)/gamma2(NS,12)
     &* xknh3(Ta)/xkelvin

c         xm12= xkNH4bb*xm23*ml6shell(I)
         xm23=xm12/xkNH4bb/ml6shell(NS)
         
         xm4= xm23+ xm12
c         if (xm4.le.1D-4) then
         xn(NS,4) =xm4/ML(1)* xn(NS,1)

                f(NS+1,4) =f(NS,4)



c     check vapour pressere
c                DO k=1,NP
c                   ml(k)= xn(i,k)*ml(1)/xn(i,1)
c                enddo
c        call vapnew(Ta,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2,pOA,pl,pl)
c        write(26,'(A,I5,5E16.6)')  'NH3', iseq,pnh3*1E6,timeeq
        
                
         endif
c         print*, 'iseqnh3, ', iseqnh3
c         print*, 'ddddd', f(ns+1,18)
cp
c     H2O equilibrium check
c     check which sells take equilibrium for varible shells for H2O
            Ishelleq=NS+1
                         
        if( imode_shell.eq.1) then
c     dlshellH2o(NSMM)
           if    (dx.le.3.0d0*dxmin.and. time.ge. .02*(x(NS+1)*1D4)**2
     & .and. dabs(RH-awshell(NS)* xkelvinH2O).le.0.01         )
     &       Ishelleq=NS
           

           NN=ns
           if ( x(NS+1)-x(NS).ge. 20E-7) Nn=NS+1
           if (imode_NH4NO3.eq.1 .or. imode_MA.eq.1) nn=1
              DO I=NS,nn,-1
              dd=dabs(f(I+1,1)-f(I,1))
c              print*,'II',i, dd,xn(I,1)
              if ( dd.le.xn(I,1) .and. I.ne.NS) then
                 Ishelleq=I+1
                 goto 291
              endif
           enddo
 291       continue


        endif

c       stop
        
         if (time.ge.timeeq) then
            Ishelleq=1
         endif


           if (imode_NH4NO3+imode_MA.ge.1 ) then
           dawmax =0.01
           timeliq = x(NS+1)**2/dh2oliq/2
           
           
           DO I=Ntr0,1,-1
c              print*, i, timetr0(I), rhtr0(I)
              
             if (timetr0(I).le.time.and. dabs( rh- rhtr0(I)).ge.dawmax)
     &              goto 236
           enddo
 236       dtimee= time- timetr0(I)
c           print*,time,rhtr0(I)
           daw=0
           DO I=1,ns
              if (dabs(RH-awshell(I)).ge. daw) daw=dabs(RH-awshell(I))
           enddo
           
           if (dtimee.ge. timeliq .and. daw.le.0.02) Ishelleq=1
        endif



        
           
         dO I=Ishelleq,ns
            f(I+1,1)= f(I,1)
         enddo

         if (idiagnose.ge.1) then
         
           print*, 'ns', ns, ishelleq
           endif
           

           
             
             
           Ta=t + fheat_11*vapheat*mm(1) * f(NS+1,1) /4/pi/rad/thermcond
           

         return
       end

      function cal_buffer(T,I,ml)
        IMPLICIT REAL*8 (A-H,m,O-Z)
       parameter (NSMM=100,np=50)
      real*8 gamma2(NSMM*2,NP)
      common /gamma2/gamma2
      real*8 ml(np)
      
      xmbuffer=0d0
        
c        cal_buffer= xmbuffer

        
c        return
c     CO2 HCO3
        xk340 = 4.448E-7*dexp(-2133*(1/T-1/298.15)) ! https://www.sciencedirect.com/science/article/pii/S0070457108703303
        xk34=xk340/gamma2(I,6)/gamma2(I,15)
        ff= xk34/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =(ml(15)+ml(13))*ff
        

c     
c     CO3-2 HCO3-
         xk2 = 10D0**(-10.33)*dexp(-3347.3*(1/t-1/298.15d0))
         xk2 = XK2*gamma2(I,15)/gamma2(I,6)/gamma2(i,14)
         ff=xk2/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(14)+ml(15))*ff


c     acetic acid



               xx= -1500.65/T-6.50923 * dlog(T)/dlog(10d0)-0.0076792* T
     & +18.67257d0
         xka= 10.d0**(xx)
         xka=xka/gamma2(I,6)/gamma2(I,8) ! takes the activity coefficient of H+
         ff=xka/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+ml(3)*ff


c     H3PO4 - H2PO4-
           xkP1= 6.9E-3*gamma2(I,24)/gamma2(I,6)/gamma2(I,25)
         ff=xkP1/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(24)+ml(25))*ff



c     HPO42- - H2PO4-

           xkp2= 6.2E-8/gamma2(I,6)/gamma2(I,26)*gamma2(I,25)
         ff=xkP2/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(26)+ml(25))*ff

c     lactic acid
        
              xka=10d0**(-3.84d0)/gamma2(I,6)/gamma2(I,33) ! take the activity coefficient of H+
         ff=xka/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(33)+ml(34))*ff
        gammah=gamma2(I,6)

c     oxalic acid
        xka=xka/gammaH/gamma2(I,30)  ! take the activity coefficient of H+
         xka2 = 5.E-5 
        xka2=xka2/gammaH/gamma2(I,31)*gamma2(I,30)  ! take the activity coefficient of OA        
         ff=xka/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(30)+ml(29))*ff
         ff=xka2/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(30)+ml(31))*ff
                       
c     NH4+, NH3
c        xkNH4bb=xknh4b(T)*gammaH/gammaNH4

        xkNH4bb=xknh4b(T)*gammaH/gamma2(I,12)
         ff=xkNH4bb*ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(4))*ff
        
c     sulfuric acid
        
        dh0=7036.629444310801
	xlk=-4.556380021818660+dh0*(1/298.15-1/T)-275./8.314*
     &  dlog(t/298.15d0)
	xks=exp(xlk)
        xks = xks /gammah/gamma2(I,28)*gamma2(I,27)
         ff=xks/ml(6)
        if (ff.ge.1) ff=1/ff
        xmbuffer =xmbuffer+(ml(28)+ml(27))*ff
        
        
        cal_buffer= xmbuffer

        return
      end
      

c       call cal_dlaw_mod(T,aw,dlacid,xvv)

      subroutine cal_dlaw_mod(T,aw,dl,xv)
       implicit real*8 (a-z)

       call cal_dlaw(T,aw,dli)
       aw0=1
           call cal_dlaw(T,aw0,dli0)
c       goto 35

c       call cal_dlaw_walker_h2o(T,aw0,dlm0)
c       call cal_dlaw_walker_h2o(T,aw,dlm)
c       eta = dlog(dlm0/dlm)/dlog(dli0/dli)

c       vmatrix=142d0
       

c       DO vv=20,2000,1D-3
c       ee= 1- 0.73*dexp(-1.79*(18d0/vv)**(1/3d0))
c       if (ee.le.eta) goto 33
c      enddoc
c 33   print*,'vv' ,vv,ee
c 35   continue
c       print*, dl, dlm,xv
       
c      if(aw.le.0.91) then
c       ff= -dlog( (1-eta) / 0.73)/1.79
c       ff=ff**3
c       vmatrix= 18/ ff
c       print*, aw,vmatrix
c      endif

      
      vmatrix=384d0
      
      
       eta = 1- 0.73*dexp(-1.79*(xv/vmatrix)**(1/3d0))
c      aw0=1d0
c      call cal_dlaw_walker_h2o(T,aw0,dlm0)
       dl = dli0* ( dli/dli0)** eta
c       print*, 'eta,dl', eta,dl
c       call cal_dlaw_walker_h2o(T,aw,dlh2o)

c       print*, dli, dl, dlH2o
       
       return
       end
      
c     radius of core in cm

      function core(time)
      implicit real*8 (a-h,o-z)
      real*8 b(10)
      data key /0/
      data keyc /0/
      logical ex
      logical ex_rcore
      real*8 timec(4000), rcore(4000), x1(1),y1(1)
      integer nc
      
      save b,ex,NB,timec, ex_rcore,rcore,nc
      
      
       INQUIRE (FILE='rcore_input.dat', exist=ex_rcore)
       if (ex_rcore) then
          if (keyc.eq.0) then
             keyc=1

             open(99,file='rcore_input.dat')
             DO I=1, 4000
                read(99,*,end=88) timec(I), rcore(i)
             enddo
 88          nc=I-1
             
             close(99)
          endif
          x1(1)=time
          N1=1
          call intpl(timec,rcore,nc,x1,y1,n1)
          core=y1(1)
          
          return
       endif
       
      
       INQUIRE (FILE='b.core', exist=ex)
       
       if (ex .and.key.eq.0) then
          key=1
          NB=0
          
          open(99, file= 'b.core')
          DO I=1,10
             read(99,*,end=10) ii, b(I)
          enddo
 10       NB=I-1
       endif

       if (ex) then
      core= b(1)              !cm
      else
         core=0d0
      endif
      
       
      
      return
      end
      

cccccccccccccccccccccccccDiffusion coefficients NH4NO3/sucrose/solution
      
       subroutine cal_dl_sucrose(T1,aw,xmv,dl)

         implicit real*8 (a-h,o-z)
         real*8 x(10)
         data key/ 0/
         save x
         if (key.eq.0) then
            CLOSE(99)
            open(99, file='b_AB_sucrose.dat')
            x(5)=0d0

                x(1)   = 0.12442E+01   
           x(2)=    0.92248E+00  
         x(3)=   -0.37361E+00  
         x(4)=     0.18181E+01  
            
            DO I=1,5
               read(99,*,end=12) II,x(I)
            enddo
 12         close(99)

            key=1
         endif
         
         t=296
         A=x(1)+x(3)*aw+x(5)*aw**2
         
         B=x(2)+x(4)*aw


!     print*,'AB', a,b
         pi=dacos(-1d0)
      xm1=18.01528
      rho1=0.999
      xv1= xm1/rho1
      rw= (xv1/4/pi*3)**(1/3d0)
      
      xm2=342.2965
      xv2= xm2/1.59
      rs= (xv2/4/pi*3)**(1/3d0)
      aw0=1d0
      call cal_dl_H2O_suc(T,aw0,dlH2o0)

            call cal_dl_H2O_suc(T,aw,dlH2o)


            
            psiw= 1- A* dexp(-B* rw/rs)
            psis= 1- A* dexp(-B)

            vma= xmv            !/6.023E23
            rma=(3./4/pi* vma)**(1/3d0)
            psima= 1- A* dexp(-B* rma/rs)


!            write(6,'(A,6E15.6)')'aw, psi',aw, psiw,psis,psis/psiw
            

                    DlMA= dlH2o0* (dlh2O/dlh2o0)**(psima/psiw)


!                    if (psis.le.0.5) Dls=1D-50
                                dlma296=dlma
!     take the temerature dependce of H2O in sucrose from Zobrist
                                call cal_dlaw_suc(T,aw,dl296)
                                call cal_dlaw_suc(T1,aw,dl)


                                dlma= dlma296* (dl/dl296)**(psima/psiw)

                                dl=dlma
                                
                                return
      
      end


c     H2O in sucrose fit Zobrist and Price data, Nadler(?)
      
      subroutine cal_dl_H2O_suc(T,aw,dl)
         implicit real*8 (a-h,o-z)

       real*8 x(3)


       x(1)=0.85094E-12
       
       x(2)=-0.12542E+01

       x(3)=0.25812E+01
c       x(3)=0.3E+01
      
       aw0=1
       call   cal_dlaw_suc(T,aw0,d0)

       
       d0=dlog(d0)
      
       
       xx=dlog(x(1))
!       xx=dlog(1D-12)
       
      
      alpha =dexp ( x(2)* (1-aw)**x(3))
      fcal = xx*(1-alpha* aw) +alpha* aw* d0
      dl = dexp(fcal)
     
      return
      end
      

c     now the binary NH3NO3 solution as a function R and aw, and xmv
      
      subroutine cal_dl_AN(T,aw,xmv,dl)

         implicit real*8 (a-h,o-z)
         real*8 x(2)
       common /xvion/ xvion
         
         pi=dacos(-1d0)
      xm1=18.01528
      rho1=0.999
      xv1= xm1/rho1
      rw= (xv1/4/pi*3)**(1/3d0)
      
      xm2=342.2965
      xv2= xm2/1.59*.999
      rs= (xv2/4/pi*3)**(1/3d0)
      aw0=1d0

            call Diff_LK(aw0,dlAN0,dlsuc)

c            psi= 1- A* dexp(-B* rw/rs)
c            psis= 1- A* dexp(-B)
c            print*,psi,psis

            call Diff_LK(aw,dlAN,dlsuc)
            dlw_dry = 4D-8 
            xx=0d0
                        call Diff_LK(xx,dlANdry,dlsuc)
c     relative to ionen
            psi =  dlog(dlandry/dlan0)/dlog(dlw_dry/dlan0)
            psi=1/psi


            vi=xvion
            vr = xmv

            rr= (vr/vi)**(1/3d0)
            A=1
            aa= (1-psi)/A

            

c            do b=0, 4,1D-6
c            ff = A*dexp(-b) -psi + 1-A*dexp(-b*rr)
!        print*, b,ff 
c            if (ff.le.0d0) goto 44
c         enddo
c 44          print*, b,ff 
            b= 0.5654
 
             DD=           A*dexp(-B) +( 1-A*dexp(-b* rr))
c            print*,'dd', rr,dd,psi
            
c            stop

            
            Dl= dlan0* (dlan/dlan0)**dd
            
                    
                    return
      
      end

cccccccccccccccccccccccccccccc
      
c     Diffusion coeffcient of sucrose in sucrose dl_suc
c     Diffusion coeffcient of NH4+ and NO3- ions in NH4NO3 solution
c     Klein 2024
      

            subroutine Diff_LK(aw,dlAN,dlsuc)
      implicit real*8 (a-h, o-z)
      a1=  5.2571978
      b1= -8.98631907
      c1=  0.75186241
      visan= a1*aw**2 + b1*aw+ c1
      visan=10d0**visan
      
      a=- 2.63720583
      b= -5.66469305
      c= 26.11320604
      d=-37.8954176
      e= 17.25426411
      vissuc= a*aw**4 + b*aw**3 + c*aw*aw + d*aw + e
      vissuc=10d0**vissuc
      t=293.15
      aw1=1d0
      vissuc1= a*aw1**4 + b*aw1**3 + c*aw1*aw1 + d*aw1 + e
      vissuc1=10d0**vissuc1
     
      
      call cal_dlaw_suc(T,aw1,dl0)
      ffs= dl0*vissuc1
      visan1= a1*aw1**2 + b1*aw1+ c1
      ffan= dl0*10d0**visan1

c      print*,'dl0',dl0,ffs
      
c      print*,'vis', vissuc,vissuc1      
      dlsuc=ffs/vissuc

      dlan=1/visan* ffan
      
      
      return
      end
      
      
c       subroutine cal_dl_sucrose(T,aw,xmv,dl)

c     diffusion coefficient in ternay NH4NO3/sucrose solution
c     ratio = molar ratio ( M_an/(M_suc+M_an)
      
      subroutine cal_dl_an_suc(T,ML,aw,xmv,dl)
      implicit real*8 (a-h,o-z)
      real*8 ML(*)
      real*8 b(5)
      save b
       data key /0/ 
      
c     treat MA as NH4NO3 at first
c     later it will be improved
      
      xman= (ml(18)+ml(12))+ml(29)
      
      ratio = (xman)/(xman+mL(2))
c      print*,'aw, ratio ', aw , ratio


      if (key.eq.0) then
         key=1
         
         b(1)=0
       b(2)=0d0
       b(3)=1d0
       b(4)=0d0
       
       open(99,file='b.mix')
       DO I=1,4
          read(99,*,end=99) ii, b(i)
       enddo
 99    continue
       close(99)
       
      endif
      
      alpha=1d0
      if (ratio.lt.0.99999d0)
     & alpha =dexp ( (b(1)+b(2)*aw+b(4)*aw**2)* (1-ratio)**b(3))

      call cal_dl_sucrose(T,aw,xmv,ds)
      call  cal_dl_AN(T,aw,xmv,dan)
      ds=dlog(ds)
      dan=dlog(dan)
       
      fcal = ds*(1-alpha* ratio) +alpha* ratio* dan
       dl =  dexp(fcal)
c       dl=0.001*dl
c       write(6,'(A,16E15.6)') 'aa ', t,aw, ds,dan,alpha,dl,xmv
c     & ,ratio       

       
      
      return
      end
      

      subroutine addbin(xn,x,NS,dd)
       implicit real*8 (a-h,m, o-z)
      parameter (NSMM=100)
      parameter (np=50, npsolid=12)         ! number of species
      real*8 xn(NSMM,np),xn0(NSMM,np),vshell(NSMM)
      common /xn/xn0


      real*8 xn2(NSMM,np),ml(np)
      real*8 x(*), x0(NSmm+1)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      
	common /m/ mm, mv,izc(NP)
      common/solid/xnsolidt

c      common /vshell/ vshell
      common /pi/ pi

       
       

       dx = x(NS+1) - x(Ns)
       dx=dx/2

       
c       if (dx.ge.dd*1.5) then

c          dx=dd
       
 
       vol = 4*pi/3d0 *( x(NS+1)**3 - (x(ns+1)-dx)**3d0)
       volshell = 4*pi/3d0 *( x(NS+1)**3 - (x(ns))**3d0)
       dd= vol/volshell
       DO J=1,NP-npsolid
          xnnn=xn(NS,j)
          xn(NS+1,j) =xnnn *dd
          xn(NS,j) =xnnn * (1-dd)
       enddo


c     nosolid split
       DO J=np-npsolid+1,np
          xn(ns+1,j)=0d0
       enddo
c
       
       x(NS+2) = x(NS+1)
       vshell(NS+1) = vol
       vshell(NS) = volshell-vol
       x(NS+1) =  x(NS+1)-dx
       NS=NS+1
       print*,ns , 'add bin'
       
c      endif
      
       return
       end     
      

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     split  when NH4+ ste > 0.3
            subroutine split_shells_NH4No3(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP),ml(np) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),xs(NSMM+1),dxliqdry(NSMM)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),vshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell,vshell

c       common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

       common /dxmin/dxmin0
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa
      real*8 m12(nsmm) ,m18(nsmm),m29(nsmm)
      x22b=0d0
      DO I=1,ns
         x22b=x22b+xn(I,2)
         vs=0
         DO J=np-npsolid+1,np
            vs=vs+ mv(j)*xn(I,j)
         enddo
         vs= vs+ 4*pi/3 * x(I)**3
         xss= ( vs/4/pi*3)**(1/3d0)
         vdry= vs+ xn(I,2)*mv(2)
         DO J=6,np-npsolid
            vdry = vdry + mv(j)*xn(I,j)
      enddo
         rdry= ( vdry/4/pi*3)**(1/3d0)
         dxliqdry(I)=rdry-xss
c         write(6,'(A,I5,4E15.6)') 'split',
c     &         i, dxliqdry(I), xn(I,12)/xn(I,1)*55.51,dxmin0
         
      enddo
c      stop
      
      dxmin=dxmin0*1.5
      
c     not possible to merge
      xmm12=0
      xmm29=0

      Do I=1,ns
         m12(I)= xn(I,12)/xn(I,1)*1000d0/mm(1)
         m18(I)= xn(I,18)/xn(I,1)*1000d0/mm(1)
         if (xmm12.le. m12(I)) xmm12=m12(I)
         m29(I)= xn(I,29)/xn(I,1)*1000d0/mm(1)
         if (xmm29.le. m29(I)) xmm29=m29(I)

      enddo
      
      if (NS.ge.50) return
      I1=0
c      n1=ns-3
c      if (n1.le.2) n1=2
      DO I=NS-1,1,-1
c         print*,'s', I, phshell(I),phshell(I+1)
         dh=dabs(m12(I)-m12(I-1))/xmm12
c         dh1=dabs(m29(I)-m29(I-1))/xmm29
c         if (dh1.ge.dh) dh=dh1
c     print*,'dh1',I, dh
         if (dh.ge..1) then

            I1=I-1
            if (vshell(I).gt. vshell(I-1)) I1=I
            I10=I1
            dx=dxliqdry(I1)
                     if (dx.le.dxmin) I1=0

             if (I1.gt.0)print*,'split NH4NO3 h1',I10,I1,dx
            
         endif

         
            dh=dabs(m12(I)-m12(I+1))/xmm12
c         print*,'dh2',I, dh
c         dh1=dabs(m29(I)-m29(I+1))/xmm29
c         if (dh1.ge.dh) dh=dh1
            
            if (dh.ge..3)  then


               I1=I+1
            if (vshell(I) .gt. vshell(I+1)) I1=I

            if (dxliqdry(I+1).ge.dxmin*2) I1=I+1
            
            dx=dxliqdry(I1)
            I10=I1
                     if (dx.le.dxmin) I1=0

             if (I1.gt.0)print*,'split NH4NO3 I+1',I10,I1,dx

            endif

            dx= dxliqdry(I1)

            if (dx.le.dxmin) I1=0
            
c            if (I1.eq.ns .and. dx.le.dxmin) I1=0
            
            if(I1.gt.0) then
            write(6,'(A,I5,5E15.6)') 'split I NH4NO3', I1,
     *m12(I1-1), m12(i1),m12(i1+1)
            xnca=0
            npl=np-npsolid
            DO II=1,NS
               print*,'before split NH4NO3', Ii,m12(Ii)
           xnca=xnca+ xn(Ii,2)
          enddo
          print*,'bfor 2', xnca

         DO I2=NS,I1+1,-1
            x(I2+2)=x(I2+1)
c     mve I1+1 to I1+2
            DO J=1,Np
               xn(I2+1,j)= xn(I2,j)
            enddo
         enddo
         
         DO J=1,NP-npsolid
            xn(I1,j)= xn(I1,j)/2
            xn(I1+1,j)= xn(I1,j)
         enddo
         DO J=NP-npsolid+1,np
            xn(I1+1,j)= 0d0
         enddo

         dx= x(I1+1)-x(I1)
         x(I1+2)= x(I1+1)
         x(I1+1)= x(I1)+dx/2
         NS=NS+1
         x22a=0d0
            DO II=1,NS
               ml(1)=1000/mm(1)
               DO J=2,np
                  ml(J)= ml(1)/xn(ii,1)*xn(ii,j)
               enddo
               IS=ii
               call calHNew(Ta,mL)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
         
          awshell(II)= aw
          phshell(II)= dlog(ML(6)*gammah)/dlog(0.1d0)
               
          print*,'after split ', Ii,awshell(Ii),xn(II,12)/xn(II,1)*ml(1)
           x22a=x22a+ xn(Ii,2)
           
        enddo
          print*,'after xn2',       x22b,x22a
          if ( dabs(x22a-x22b)/x22a.ge.1D-6) then
             print*, ' split mass not conserved stop'
             stop
          endif
          

         return
         endif

      enddo
      
c      print*,'NS split',ns
      

          return
          end

      
      subroutine merge_NH4NO3(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np),xna(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),ml(np),dxshell(nsmm)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),m12(NSMM),m29(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell
     & ,vshell0(nsmm)
c     common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa
 
      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3,
     &       imode_MA
 
      common /dxmin/dxmin0,factorshell      
      common /NS0/NS0
 22   continue
      
c     merge or adjust top layer
      if (NS.le.1) return
c      print*,'meger NH4NO3'
      
      vv=x(1)**3
      x22b=0d0
      DO I=1,ns
         vv1= xn(I,1)*mv(1)+xn(I,2)*mv(2)
         DO J=7, np! -npsolid
            vv1=vv1+ mv(j)*xn(i,j)
         enddo
         vv1= vv1/4/pi*3d0  !r**3
         dxshell(I) = (vv + vv1) **(1/3d0) - vv**(1/3d0)
         vv=vv+vv1
         m12(I)= xn(I,12)*55.51/xn(i,1)
         m29(I)= xn(I,29)*55.51/xn(i,1)
         x22b=x22b+xn(i,2)
c                  print*,'mm',i, dxshell(i), m12(I)
      enddo

c     check maximing thickness jump a facttor of 5
      factor = factorshell
      
 
      
      ismerge=0

            if (dxshell(NS-1).le. dxshell(NS)/3) then

            print*,'merger NS-1',NS,dxshell(NS-1), dxshell(NS)
               
               DO J=1,Np-npsolid
                  dd=(xn(ns-1,j)+xn(ns,j))/2d0
                  xn(NS-1,j)= dd
                  xn(NS,j)= dd
               enddo
               
               return
            endif


            
      DO I=1,NS
         dd=dxmin0*.5
c         if (I.eq.NS) dd=.2D-7
         if (I.lt.NS)dd=dxmin0*.75
         if (dxshell(I).le.dd ) then

c
c            merge or increase the thickness
c     merge with the shell below
            print*,'aaaa', I,NS
            
            if (I.gt.1) then
c               if ( dxshell(I-1)+dxshe 

               if  (dxshell(I)+dxshell(I-1).le.2.1*dxmin0) then


                  write(6,'(A,I5,6E15.6)') 'Merge NH4NO3 down' , I,x22b
                  DO J=1,Np
                     xn(I-1,j)=xn(I-1,j)+xn(I,j)
                  enddo
                  
                  DO II=I,NS-1
                  DO J=1,Np
                     xn(Ii,j)=xn(Ii+1,j)
                  enddo
                     enddo

                    NS=NS-1
                     x22a=0d0
                     DO II=1,ns
                        x22a=x22a+xn(Ii,2)
                     enddo
                     write(6,'(A,I5,6E15.6)') 'after' , I,x22a
                     
               else
                  write(6,'(A,I5,6E15.6)') 'Merge NH4NO3 rebin' , I

                  ff=0.5d0
c                  if (NS.eq.I)then
c                  vmin= x(NS+1)**3-(x(ns+1)-dxmin0)**3
c                  ff=vmin*1.2/( x(NS+1)**3-x(ns-1)**3)
c                  if (ff.ge.0.5d0) ff=0.5                  
c               endif

                  
               print*,'ff', ff
                  DO J=1,np-npsolid
                        xx=xn(I,j)+xn(I-1,j)
                        xn(i,j)=xx*ff
                        xn(i-1,j)=xx*(1d0-ff)
                     enddo
                     
c                     stop
                     
                  endif
            ismerge=1
                  goto 33
            endif

c c            stop
            
            if (I.lt.NS) then
c               if ( dxshell(I-1)+dxshe 

               if  (dxshell(I)+dxshell(I+1).le.2.1*dxmin0) then
                  write(6,'(A,I5,6E15.6)') 'Merge NH4NO3 up' , I

                  DO J=1,Np
                     xn(I,j)=xn(I,j)+xn(I+1,j)
                  enddo
                  
                  DO II=I+1,NS
                  DO J=1,Np
                     xn(Ii,j)=xn(Ii+1,j)
                  enddo

                     enddo

 331                 NS=NS-1
         goto 33            
            ismerge=1
                  endif
               
            endif
            
          
         endif

c     merge up

         
         
      enddo
      
 33   continue
c      print*,'NS', NS
      
c      read(5,'(A1)')
c     merge NH4NO3
      if (NS.le.NS0+1 .and. awshell(NS).le.0.2) return
      if (NS.le.NS0 .and. awshell(NS).ge.0.2) return
c     Merger 1, NS-2
      DO I= 1, Ns-4
         ism=1
         if (dxshell(I).ge. dxshell(I+1)*factor*0.5 ) ism=0
         if (dxshell(I+1).ge. dxshell(I+2)*factor*.4 ) ism=0



         if (dabs(awshell(I)-awshell(I+1)).ge.0.003 ) ism=0
         if (imode_NH4NO3.eq.1  ) then
            if ( dabs( m12(I)-m12(I+1))/m12(I) .gt. .02 )ism=0
         endif
         if (imode_MA.eq.1) then
            if ( dabs (m29(I)-m29(I+1))/m29(I) .gt. .02 )ism=0
         endif
         if (imode_pH.le.1) then
            if ( dabs( phshell(I)-phshell(I+1)) .gt. .02 )ism=0
         endif
      
      
         if (ism.eq.1) then
            DO J=1,Np
               xn(I,j) = xn(I,j)+ xn(I+1,j)
            enddo
               DO II=I+1,NS-1
               DO J=1,Np
                xn(Ii,j) =  xn(Ii+1,j)
               enddo
                  enddo

         print*,'ddd', dxshell(I+1), dxshell(I+2)*factor*.6 

         print*, 'dd ',dxshell(I),dxshell(I+1),dxshell(I+2)
                  NS=NS-1
            
            print*, 'mergen I,I+1 aw,NS', I,I+1,NS
            ismerge=1
            
            goto 333
         endif
         

      enddo
 333  continue
      if (ismerge.eq.1) then
      x22a=0d0
      DO I=1,NS
         x22a=x22a+xn(I,2)
      enddo
      print*,' merge xn22', x22b, x22a
      if ( dabs(x22a-x22b)/x22a.ge.0.01) then
         print*,'merge, mass not conserved , stop'
         stop
      endif
      endif


      if( ns.lt.NS0) then
c     split the shell with the highst volume
         xx=0
         DO I=1,NS
            if (xn(I,1).ge.xx) then
               ii=I
               xx=xn(I,1)
            endif
         enddo
c     split II
         DO I=NS,II+1,-1

            DO J=1,np
               xn(I+1,j)=xn(I,j)
            enddo
         enddo
         
            DO J=1,np
               xx= xn(II,j)
               xn(II,j)=xx/2
               xn(II+1,j)=xx/2
            enddo
            
            
            NS=nS+1
            
         endif
         
      

      
      return
      
      

          end



c split factor
     
      subroutine split_factor(x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np),xna(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),ml(np),dxshell(nsmm)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),m12(NSMM),m29(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell
     & ,vshell0(nsmm)
c     common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa
 
      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3,
     &       imode_MA
 
      common /dxmin/dxmin0,factorshell      
      common /NS0/NS0
 22   continue
      
c     merge or adjust top layer
      if (NS.le.1) return
c      print*,'meger NH4NO3'
      
      vv=x(1)**3
      x22b=0d0
      DO I=1,ns
         vv1= xn(I,1)*mv(1)+xn(I,2)*mv(2)
         DO J=7, np! -npsolid
            vv1=vv1+ mv(j)*xn(i,j)
         enddo
         vv1= vv1/4/pi*3d0  !r**3
         dxshell(I) = (vv + vv1) **(1/3d0) - vv**(1/3d0)
         vv=vv+vv1
         m12(I)= xn(I,12)*55.51/xn(i,1)
         m29(I)= xn(I,29)*55.51/xn(i,1)
         x22b=x22b+xn(i,2)
c                  print*,'mm',i, dxshell(i), m12(I)
      enddo

c     check maximing thickness jump a facttor of 5
      factor = factorshell
      
      nn=1
      DO I=Ns,nn,-1
         dd=dxshell(I)
         if (dd.le.dxmin0) dd= dxmin0
         if  (dxshell(I-1) .ge. dd*factorshell) then
c     split bin I-1
            print* , 'Split jump', I-1,NS

            DO II=NS,I,-1
               DO J=1,NP
                  xn(II+1,j)=xn(II,j)
               enddo
            enddo
            
               DO J=1,NP
                  xn(I-1,j)=xn(I-1,j)/2
                  xn(I,j)=xn(I-1,j)

               enddo

               NS=NS+1
               
            goto 22

            
            
         endif
         
      enddo
      
 
      return
      
      

          end

      

c     reset_shells_NH4
      subroutine reset_shells_NH4NO3(NS,x)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=50)
      parameter (NSMM=100,npsolid=12)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np),xna(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM),ml(np),dxshell(nsmm)
      
      common /m/ mm, mv,izc(NP)
      common /pi/pi
      common /xn/xn
      real*8 ml6shell(NSMM),awshell(NSMM),sshell(NSMM,npsolid)
      real*8 xhshell(NSMM),phshell(NSMM),m12(NSMM),m29(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sshell,phshell
     & ,vshell0(nsmm)
c     common/flux/T,Ta
      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl, partvapOA,partvaplac

      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa

      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3,
     &       imode_MA
 
       common /dxmin/dxmin0
      common /NS0/NS0

c     merge or adjust top layer
c     if (NS.le.ns0) return
c      print*,'meger NH4NO3'
      is=1
      
      DO I=1,NS-1
         if (dabs( awshell(I)-awshell(I+1)).ge.0.005) is=0d0
      enddo
c      print*,'is', is
c      DO I=1,NS
c         xn(I,2)
c        enddo
      DO I=1,NS-1
         ff=dabs( xn(I,2)/xn(I,1)-xn(I+1,2)/xn(I+1,1))/
     &        (xn(I,2)/xn(I,1))
         dd=0.1
         if (I.ge.NS-3) dd=.03
         if (dabs(ff).ge.dd) is=0
c         print*,i, xn(I,2)/xn(I,1)*55.51,ff
         
      enddo
c      print*,'is', is

      DO I=1,NS-1
         ff=dabs( xn(I,12)/xn(I,1)-xn(I+1,12)/xn(I+1,1))/
     & (xn(I,12)/xn(I,1))    
         dd=0.1
         if (I.ge.NS-3) dd=.03

         if (dabs(ff).ge.dd) is=0

      enddo
      
      if (is.eq.1) call reset_shells(NS,x)
c      print*,'is', is
      

      
      return
      
          end

      


	    function RI_ANSuc(wt_suc,wt_AN,xl,T,rhoo)
	    implicit real*8 (a-h,o-z)
        real*8 a_mix, x, wt_suc, wt_AN, wt 
        real*8 x_h2o, x_suc, x_AN, test
        real*8 mol_h2o, mol_suc, mol_AN, total_mol
C23456	called funktions: AH2O, a_AN, a_suc

c        conversion from weight fraction to mole fraction         
         wt = wt_suc+ wt_AN
         mol_h2o = (1-wt)/18.016D0
         mol_suc = wt_suc/342.33d0
         mol_AN = wt_AN/(18d0+62.01d0)
         
         total_mol =  mol_h2o + mol_suc + mol_AN
    
         x_h2o = mol_h2o/total_mol
         x_suc = mol_suc/total_mol
         x_AN = mol_AN/total_mol
         
         test = x_AN + x_h2o + x_suc

        x = x_suc + x_AN
        A = (1-x)*ah2O(xl)+ x_suc*a_suc(wt,xl,T,rhoo)
     &        + x_AN*a_AN(wt,xl,rhoo)

c        print*,a_suc(wt,xl,T,rhoo),a_AN(wt,xl,rhoo)
c        print*, 'aa',x,ah2o(xl)
        
        a_mix = A * rhoo/
     *   ((1-x)*18.016D0+
     *   x_suc*342.33d0+
     *   x_AN*(18d0+62.01d0))       
        RI_ANSuc=dsqrt((1.D0+2.D0*a_mix)/(1-a_mix))
c        print*,'a', x_suc, x_an,rhoo,a,a_mix
        
        
c        print *, a_AN(wt,xl,rhoo), ah2O(xl), RI_ANSuc
c        print*,' wt_an,wt_suc', wt_an, wt_suc

        
        return
	end
c-------------------------------------------------------------------
      
       function rhomm(m,mm,mv,ns,xm)
       implicit real*8 (a-h,k,m, o-z)
       real*8 m(*),mm(*), mv(*)
       xm=0
       xv=0

       DO J=1,NS
          xm=xm+m(J)*mm(J)
          xv=xv+m(J)*mv(J)
                 enddo
                 rhomm= xm/XV
                                  return
       end

	function a_suc(wt, xl, TK,rhoo)  
c     wt in fraction, T in Kelvin, xl in um; --> n in air
        implicit real*8 (a-h,o-z)
        real*8 x(6)
        real*8 n, Pii, tau, xl, Lambda, nair, n_suc, a, aw, wt
        real*8 x_suc, x_h2o, AA, tk, t
        real*8 mol_suc, mol_h2o, total_mol, test
        parameter (NP=50) 
        real*8 mm(NP),mv(np),m(np)
        
         common /mm/ mm,mv
        
        nair=1.000272
        Pii=(wt*100d0-35)/35
        tau=(TK-273.15-28)/10
        Lambda=((xl*1000)-546)/43

c23456
       n=1.3911588-0.00157931*tau+0.065314306*Pii
     *       -0.0000931172*tau**2+0.0090056*Pii**2
     *       -0.000438616*Pii*tau+0.000841194*Pii**3
     *       +0.0000523892*Pii*tau**2
     *       -0.00005514751*Pii**4+0.0000262176*Pii**2*tau
     *       -0.0000471533*Pii**5+0.0000377621*Pii**3*tau
     *       -0.00165484*Lambda+0.0000098966*Lambda*tau
     *       -0.000209707*Lambda*Pii
     *       -0.00002775*Lambda*Pii**2

         n_suc = n/nair
         
c         wt = wt/100
         a  = (n_suc**2-1)/(n_suc**2+2)
          
         mol_suc = wt/342.33d0
         mol_h2o = (1-wt)/18.016D0
         total_mol = mol_suc + mol_h2o
    
         x_suc = mol_suc/total_mol
         x_h2o = mol_h2o/total_mol
c        test = x_suc + x_h2o
c        print*, wt, mol_h2o, total_mol, x_h2o, test
c     23456
c         m=0d0
c         m(1)= 1000d0/mm(1)
c         m(2) = m(1)/x_H2O* x_suc
c         rhoo= rho_mix(tk,m,mm,mv)

                  AA = (x_suc*342.33d0+x_h2o*18.016D0)/rhoo*a
                  
c     *  /rhosuc(wt,tk)*a
         
       a_suc = (AA - x_h2o*ah2o(xl))/x_suc ! Molar refractivity of pure sucrose 
c       wt100=wt*100
c       write(6,'(A, 3F12.4)')
c     &    'rho, a_suc,rhosuc', rhoo, a_suc ,rhosuc(wt100,tk)
        

c     print*, AA, x_h2o, x_suc, a_suc, ah2o(xl)
c       print*, AA, a_suc, wt
        return 
        end



c-------------------------------------------------------------------
C	Molar refractivity of NH4NO3 (from Tang 1981) probably sodium D line (589 nm)
c-------------------------------------------------------------------
      function a_AN(wt, xl,rhoo) 
        implicit real*8 (a-h,o-z)
        real*8 wt, n, a,t, aw, AA, mol_AN, mol_h2o,total_mol, x_AN, x_h2o
        real*8 x(4)
        parameter (NP=50)
        real*8 mm(NP),mv(NP),m(NP)
        common /mm/ mm,mv
        
      data x
     */    1.3337,  0.119,  1.3285, 0.145/

      if(wt .le. 0.205) then
       n= x(1) + x(2)*wt
      else
       n = x(3) + x(4)*wt
      endif
      
      
       a = (n**2-1)/(n**2+2)

       
         mol_AN = wt/(18d0+62.01d0)
         mol_h2o = (1-wt)/18.016D0
         total_mol = mol_AN + mol_h2o
    
         x_AN = mol_AN/total_mol
         x_h2o = mol_h2o/total_mol
c        test = x_AN + x_h2o
c23456   print*, wt, mol_h2o, total_mol, x_h2o, test
c         m=0d0
c         m(1)= 1000d0/mm(1)
c         m(3) = m(1)/x_H2O* x_AN
c         m(4)=m(3)
c         tk=298.15d0
c         rhoo= rho_mix(tk,m,mm,mv)
c         print*,'rhoo', rhoo
         
         AA = (x_AN*(18d0+62.01d0)+((1-x_AN)*18.016D0))/
     *    rhoo*a ! Molar refractivity of the binary Ammoniumnitrate water solution 

         a_AN = (AA - x_h2o*ah2o(xl))/x_AN ! Molar refractivity of pure Ammoniumnitrate 
c         print*, 'rho_an, a ', rhoo, a_an
         
c        write(18,*) wt, AA, a_AN
        return
        end

      function rho_mix(t,m,mm,mv)
       implicit real*8 (a-h,k,m, o-z)
       parameter (np=50)
              real*8 m(np),mm(np), mv(NP)

              xm=0
       xv=0

c       call vapnew(T,M,aw,pnh3,pHNO3,gammaH,gammaNO3,gammaNH4)
c       call cal_mv(aw,t,mv)

       
       DO J=1,Np
          xm=xm+m(J)*mm(J)
          xv=xv+m(J)*mv(J)
        enddo
                 rho_mix= xm/XV
                                  return
       end

	function ah2o(lambda)
	implicit real*8 (a-h,o-z)
	real*8 b1(5),b2(5),lambda
	data b1 /1.4211,7.540441,-9.81512155,5.7218513,-1.237611/
c	data b2/3.506813,.16411,-.0461273,.012763141,-3.579466E-4/

       data b2
     * /3.5740152,5.820734E-2,1.35319E-2,-1.5603E-3,8.85744E-4/
     


	x=1./lambda
	if(x.le.1.23) then
	ah2o=b1(1)+b1(2)*x+b1(3)*x**2+b1(4)*x**3+b1(5)*x**4
	else
	ah2o=b2(1)+b2(2)*x+b2(3)*x**2+b2(4)*x**3+b2(5)*x**4
	endif
	return
	end


      

      subroutine cal_h2o_guess(T,aw,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1), x(3)

      data key /0/
      save x
      
      if (key.eq.0) then

      key=1
      

         open(99, file='b_H2O.dat')

         x(1)=  0.14021E-08
         x(2)=-0.19276E+01 
         x(3)= 0.69964E+00
         
         read(99,*,end=23)  xx,x(1)
         read(99,*,end=23) xx,x(2)
         read(99,*,end=23)xx, x(3)

         

 23      continue
         close(99)

      endif
      
        aw0=1
       call   cal_dlaw_citric(T,aw0,d0)
       t298=298.15
       call   cal_dlaw_citric(T298,aw0,d298)
       d0=d0*2.44E-5/d298
       
       
       
       d0=dlog(d0)
      
       
       xx=dlog(x(1))
!       xx=dlog(1D-12)
       aw0=0d0
       call cal_dlaw_citric(T,aw0,dlt)
       aw0=0d0
       t0=298.15
       call cal_dlaw_citric(T0,aw0,dlt0)

       xx=xx+ dlog( dlt/dlt0)
       
       alpha=1
       if (aw.lt.1d0)       alpha =dexp ( x(2)* (1-aw)**x(3))
      fcal = xx*(1-alpha* aw) +alpha* aw* d0
      dl = dexp(fcal)

      
       return
      end
      


      subroutine cal_dlaw_eta(T,aw,dl,xv)
       implicit real*8 (a-z)
       real*8 x(5)
       data key /0/
       save x
       if (key.eq.0) then
          key=1
          CLOSE(99)
          open(99,file='b_AB_SLF.dat')


          x(1) = 0.23592E+01 
          x(2)=-0.10941E+01 
          x(3)=3d0 
          x(4)=0.54564E+00    
          x(5)=150d0
          DO I=1,5
             read(99,*,end=22) xx, x(I)
             enddo
 22          continue
             
          endif
       
       call cal_h2o_guess(T,aw,dlH2o)

       aw0=1
       call cal_h2o_guess(T,aw0,dlH2o0)
c       print*, t, aw0, dlh2o0
       
c       stop
       vmatrix=x(5)
       

       A=x(1)+ x(3)*aw
       B=x(2)+ x(4)*aw
       
         etaw = 1- A*dexp(-B*(18d0/vmatrix)**(1/3d0))
         eta = 1- A*dexp(-B*(xv/vmatrix)**(1/3d0))


       dl = dlh2o0* ( dlh2o/dlh2o0)** (eta/etaw)


       
       return
       end
      
      
            function aw_corr(aw)
      implicit real*8 (a-h,o-z)
      real*8 x4(5),y4(5) ,x1(1),y1(1)
      data x4 / 0d0, .6d0, 0.75, .85d0, .99d0/
      data y4 / .87d0, .87d0, 1d0,1.06d0, 1d0/
      integer N1,n4

      N1=1
      N4=5
      x1(1)=aw
      call intpl(x4,y4,n4,x1,y1,n1)
c      write(6,*) x1(1),y1(1)

      aw_corr=aw* y1(1)
      
      
      return
      end

      subroutine cal_MV(aw,t,mv)
       implicit real*8 (a-h,m, o-z)
       parameter (np=50)
       real*8 mv(np)

      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
     & , imode_MA,imode_EDB
      

       mv(1) = f_v1(T)
       if (imode_NH4NO3.eq.1) then
          mv(2)=mv_suc(aw,t)
       endif
       
       mv(12)=mv_NH4(aw,t)
       mv(18)=mv_NO3(aw,t)
        return
       end


      function mv_NH4(aw,t)
        implicit real*8 (a-h,m,o-z)
         mv_NH4 = mv_NH4NO3(aw,t)- mv_NO3(aw,t)
        return
      end
 
      function mv_NH4NO3(aw,t)
        implicit real*8 (a-h,m,o-z)

        real*8 x(14)
        data x/
     1   -0.88414E+02,
     2    0.52078E+00,
     3   -0.16461E-03,
     4   -0.13216E+00,
     5   -0.26002E-01,
     6    0.91500E-04,
     7    0.10380E-01,
     8   -0.27435E-01,
     9   -0.29104E-04,
     1   -0.65107E+01,
     1    0.37660E+00,
     1   -0.92953E-03,
     1   -0.10092E+03,
     1    0.28167E+00 /


c     density data of HNO3 https://www.handymath.com/cgi-bin/nitrictble2.cgi?submit=Entry        
      c0= x(1)+ x(2)*t+ x(3)*t**2
      c1= x(4)+ x(5)*t+ x(6)*t**2
      c2= x(7)+ x(8)*t+ x(9)*t**2
      c3= x(10)+ x(11)*t+ x(12)*t**2
      c4= x(13)+ x(14)*t


      mv_NH4NO3= c0 + c1*aw + c2*aw**2 + c3*aw**3+ c4*aw**4

        return
      end
      
      
      function mv_NO3(aw,t)
        implicit real*8 (a-h,m,o-z)
        real*8 x(12)
        data x/
     1    0.10947E+02,
     2    0.12881E+00,
     3   -0.10524E-03,
     4   -0.18117E+01,
     5   -0.25917E+00,
     6    0.61919E-03,
     7    0.19963E+02,
     8    0.24979E+00,
     9   -0.79861E-03,
     1   -0.97904E+02,
     1    0.44167E+00,
     1   -0.49184E-03/

c     density data of HNO3 https://www.handymath.com/cgi-bin/nitrictble2.cgi?submit=Entry        
      c0= x(1)+ x(2)*t+ x(3)*t**2
      c1= x(4)+ x(5)*t+ x(6)*t**2
      c2= x(7)+ x(8)*t+ x(9)*t**2
      c3= x(10)+ x(11)*t+ x(12)*t**2

      xvHplus=1d0               ! molar volume H+
      mv_NO3= c0 + c1*aw + c2*aw**2 + c3*aw**3-xvHplus
      return
      end

            function mv_suc(aw,t)
        implicit real*8 (a-h,m,o-z)
        real*8 x(12)
        data x/
     1    0.14686E+03 ,
     2    0.48408E+00 ,
     3   -0.78302E-03 ,
     4   -0.14731E+01 ,
     5   -0.48379E+00 ,
     6    0.14771E-02 ,
     7   -0.88865E+00  ,
     8    0.11680E+01  ,
     9   -0.35647E-02  ,
     1   -0.99821E+02  ,
     1   -0.21440E+00  ,
     1    0.15593E-02  /      
        
      c0= x(1)+ x(2)*t+ x(3)*t**2
      c1= x(4)+ x(5)*t+ x(6)*t**2
      c2= x(7)+ x(8)*t+ x(9)*t**2
      c3= x(10)+ x(11)*t+ x(12)*t**2
      mv_suc= c0 + c1*aw + c2*aw**2 + c3*aw**3

      return
      end
      
      

            function f_v1(t)
      implicit real*8 (a-z)
      tt=t-273.15d0
      rh0=0.99989d0+ 0.52785E-04*tt-0.72916E-05*tt**2+tt**3*0.31399E-07  
      xm1=18.01528 
      f_v1= xm1/rh0*18.0152797698975d0/18.0463148620371d0
      

      return
      end
      

      subroutine aw_back
     & (t,MLa,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

      implicit real*8 (a-h,m,o-z)
      integer NP
      parameter (np=50, npsolid=12)
      real*8 ML(Np),mla(*),ml0(NP)


      real*8 MM(NP) ! molar mass
      real*8 Mv(NP) ! molar volume
      integer izc(NP)

      common /M/ MM,mv          ,izc
      common /awin/ awin,aws,xvol,xmi
      data key /-1/
      
c      parameter (NM2=20002)
c      real*8 aw22(nm2,nm2) ,gno32(nm2,nm2),gnH42(nm2,nm2)

      parameter (NM2=50002,nt=101)

      real*8 gNH4(nt,NM2)
      real*8 gNO3(nt,NM2)
      real*8 awin2(nt,NM2)
      real*8 aws2(nt,NM2)
      

      common/solid/xnsolid
c      common /output/ imode_output,idiff,imode_pH,imode_eq
      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
     & , imode_MA      
      

      real*8 gh(nt,NM2)
      real*8 gOA(nt,NM2)
      real*8 gHOA(nt,NM2)
      
      save data,key, awin2,aws2,gNH3,gNH4,gh, goa, ghoa
        common/gammas/ gammas1,gammas2, gammaHCO3, gammaco3,gammah2po4
     &,gammahpo4 ,gammaHOA,gammaOA,gammah3po4,gammaca,gammamg,gammak
     &     , gammah2oA
      

       
       n2=2
       n1=1
c     no look up take for imode_output 4

       xms=mla(28)+mla(27)+mla(29)+mla(30)+mla(31)

c       print*, 'imode ', imode_NH4NO3
       
c       if(imode_output.eq.4.or. xms.gt.1d-30.or. imode_output.eq.3) then
       if (imode_NH4NO3 + imode_MA .eq.0) then
       call aw_back_model
     & (t,MLa,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      else
         if (key.le.-1) then
            key=1
c     p
            print*,'generatinf lookup table'

            tt=288
            I=0
            DO Itt=1,nt
               tt= 280+ .2*(Itt-1)
               print*,ITT
            if (imode_NH4NO3.eq.1) then   
            DO I=1,NM2
               ML=0
            ml(1)=1000d0/mm(1)
               xmm=1D-5 + (I-1)/200d0
                  ml(2)=0d0
               ml(4)= xmm
               ml(12)= xmm
               ml(18)= xmm
      call aw_back_model
     & (tt,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
               awin2(itt,I)= aw
               gNO3(itt,I)= gammaNO3
               gNH4(itt,I)= gammaNH4
               ML=0
            ml(1)=1000d0/mm(1)
            ML(2)=xmm
      call aw_back_model
     & (tt,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
               aws2(itt,I)= aw
               

            enddo
         endif

c     aicd
            if (imode_MA.eq.1) then   
               nm22=5002
               DO I=1,NM22
               ML=0
            ml(1)=1000d0/mm(1)
            xmm=1D-5 + (I-1)/20.
               if (itt.eq.1) then
c            print*,'xmm',i,xmm
            endif
            ml(2)=0d0
               ml(29)= xmm

      call calhnew_model(tt,ml)
      call aw_back_model
     & (tt,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
               awin2(itt,I)= aw
               gh(itt,I)= gammah
               ghoa(itt,I)= gammahoa
               goa(itt,I)= gammaoa
               ML=0
            ml(1)=1000d0/mm(1)
            xmm=1D-5 + (I-1)/20.
            ML(2)=xmm
      call aw_back_model
     & (tt,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
               aws2(itt,I)= aw
               if (itt.eq.1) then
                  write(6,'(2I5,15E15.6)')itt,i,
     &                  xmm,ml(29),ml(30),ml(31),gh(itt,i)
     &              , ghoa(itt,i), goa(itt,i),aws2(Itt,i),awin2(Itt,i)
                    
               endif
               

            enddo
         endif
         
         enddo
            print*, I,J

 
            
            print*,'finish lookup table'
c         stop
            
         endif
         Itt= (t-280)/0.2d0+1.5d0
         if (itt.le.1) itt=1
         if (itt.gt.nt) itt=nt
           gammah=1d0
           gammas1=1d0
           gammas2=1d0
           gammaNH4=1d0
           gammaNa=1d0
           gammacl=1d0
           gammaNO3=1d0
           gammaH3pO4=1d0
           gammaHOA=1d0
           gammaOA=1d0
           gammaCO3=1d0
           gammaHCO3=1d0
           gammaK=1
           gammaca=1
           gammamg=1
           dx=1/200d0
           if (imode_MA.eq.1) dx=1/20d0
           
            I= (mla(2)-1D-5)/dx +1.5d0
            if (I.ge.NM2) I=nm2
            if (I.le.1) I=1
            aws=aws2(itt,I)

             if (imode_nh4no3.eq.1) then
             dx= 1/200d0

            dd=(mla(12)+mla(18))/2d0
            
            I= (dd-1D-5)/dx +1.5d0
            
            if (I.ge.NM2) I=nm2
            if (I.le.1) I=1
            awin=awin2(itt,I)
            gammaNO3= gno3(itt,I)
            gammaNH4= gnH4(itt,I)
         endif
         
c
           if (imode_MA.eq.1) then
             dx= 1/20d0

             dd=(mla(29)+mla(30)+ml(31))
  
            I= (dd-1D-5)/dx +1.5d0

            nm22=5002
            if (I.ge.NM22) I=nm22
            if (I.le.1) I=1
            awin=awin2(itt,I)
            gammah= gh(itt,I)
            gammaHOA= ghoa(itt,I)
            gammaOA= goa(itt,I)
         endif
         
c 
c          print*, i, gammaNH4
c            if (I.ge.1000)stop
            aw=aws*awin
            
c            print*, aws,awin,aw
            


      endif
      if (imode_EDB.ne.1) return

       if (imode_EDB.eq.1) then
        aww = aw
        aw1=.75
        aw2=.55
        
        if (aww.ge.aw1) aww = aw1
        if (aww.le.aw2) aww = aw2
                ff = 1 -.11 * ( aw1- aww )/ (aw1-aw2)

                aw = aw*ff
                
      endif

      return

      end
      


	function gammasn(T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

	implicit real*8 (a-h,o-z)

	Parameter( nmax=12 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax),c1(nmax,nmax)
	real*8 Bphi(nmax,nmax),C(nmax,nmax)
	real*8 I, I2, MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
        real*8 B20(nmax,nmax),alpha2(nmax,nmax)
	common /alpha/ alpha2, b20
        common /xm24/ xm24,xlamc,xlama1,xlama2,xlamnn
        common /xm29/ xm29,xlam1_29,xme1_29,xlam2_29,xme2_29
       logical ex
       common /ex/ ex


       fphi1=0d0
       xs=0d0
       f3=0d0
       f4=0d0
       f5=0d0
       f6=0d0              
  
	

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
           if (mC(IC).gt.0) then
           I=I+mC(IC)*zC(IC)**2
           endif
	enddo
	DO IA=1,NA
           if (mA(IA).gt.0) then
	I=I+mA(IA)*zA(IA)**2
        endif
      enddo
	I=I/2.
        if (I.le.0d0) then
           gammasn=1d0
           goto 333
        endif

	I2=sqrt(I)

CCC	Calculate Bphi, C
	x=sqrt(I)*alpha

	DO Ic = 1, NC
	DO Ia = 1, Na
           
	   if (alpha2(ic,ia).le.0d0) alpha2(ic,ia)=1

           if (MA(IA).gt.0d0 .and. Mc(Ic).gt.0d0) then

      	   x2=sqrt(I)*alpha2(ic,ia)
	   Bphi(IC,IA) = b0(Ic,Ia) +exp(-x)*b1(Ic,Ia)+exp(-x2)*b20(Ic,Ia)
c           Bphi(IC,IA) = b0(Ic,Ia) +exp(-x)*b1(Ic,Ia)
        omega1=omega(Ic,IA)
        xo=omega1*I2
	C(Ic,Ia) = C0(Ic,Ia)+c1(Ic,IA)*exp(-xo)
        endif
      enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2
C	calculate sum mi and Z

	z=0	
	xmi=0
	DO IC=1,NC
           if (mC(IC).gt.0) then           
           Z=Z+ZC(IC)*MC(IC)
           xmi=xmi+MC(IC)
           endif
	enddo
	DO IA=1,NA
           if (ma(ia).gt.0) then

           Z=Z+ZA(IA)*MA(IA)
	xmi=xmi+MA(IA)
        endif
      enddo


	fphi1=  -Aphi*I**(3./2.)/(1+1.2*I2)

	xs=0
	DO Ic=1,NC
	DO Ia=1,Na
           if (MA(IA).gt.0d0 .and. Mc(Ic).gt.0d0) then
           xs=xs+MA(IA)*MC(IC)*(Z*C(IC,IA)+Bphi(IC,IA))
           endif
	enddo
	enddo

	F3 = 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
c	z1=ZC(IC1)
c     z2=ZC(IC2)
           if (MC(Ic1).gt.0d0 .and. MC(IC2).gt.0d0) then
	J1=ZC(IC1)+.1
	J2=ZC(IC2)+.1
	IF(j1.eq.j2) goto 21
	call EFUNC(J1,J2,Aphi,I,E,ED)
	pp=e+I*ed	

	F3=F3+pp*MC(Ic1)*MC(IC2)
21	continue
        endif
	enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
c	z1=ZA(IA1)
c	z2=ZA(IA2)
           if (MA(IA1).gt.0d0 .and. MA(IA2).gt.0d0) then
	J1=ZA(IA1)+.1
	J2=ZA(IA2)+.1
	IF(j1.eq.j2) goto 22
	call EFUNC(J1,J2,Aphi,I,E,ED)
	pp=e+I*ed	
	F4=F4+pp*MA(IA1)*MA(IA2)
22	continue
        endif

      enddo
	enddo
	xmh0=0
	DO IA=3,NA
	xmh0=xmh0+mA(IA)
	enddo

        f5 =xm24*xm24 /2 * xlamnn + mc(1)*xlamc*xm24 + ma(7)*xm24*xlama1
     &       + ma(8)* xm24*xlama2

         aa3= dsqrt(xm29/2)
       xlam_29=xlam1_29*dexp(-aa3/xmE1_29)+xlam2_29*dexp(-aa3/xme2_29)
       f6 =     xm29*xm29 /2 * xlam_29
c       print*, 'xm29', xm29
c       print*,'xlam1', xlam1_29
c       print*,'xme1', xme1_29

c       print*,'xlam2', xlam2_29
c       print*,'xme2', xme2_29
c       stop
       
 333   continue
       if (ex) then
          
        phix=  fphi1 + xs  +f3+f4+f5+f6
	phi=1+ phix*2/(xmi+xm24+xm29)
	as=(-phi)*18d0/1000.*(xmi+xm24+xm29)
        gammasn=exp(as)
      else

         phix=  fphi1 + xs  +f3+f4+f5
	phi=1+ phix*2/(xmi+xm29)
	as=(-phi)*18d0/1000.*(xmi+xm29)
        gammasn=exp(as)
         endif
         
	return
	end
      

      subroutine cal_rcore(time,  rliq,rcore)
      
      
      implicit real*8 (a-h,o-z)
      
      real*8 timeIM(3000),v(3000),rIM(3000)

      character*20 ti

      data key /0/
      real*8 x1(1),y1(1)
      logical ex_IM
      
      
      save TImeIM, NIM, rIM,ex_IM,key

      common /fvoll/ rm,rliqq
      external fvol
      rliqq=rliq
      
      
      if (key.eq.0 ) then
         key=1
         INQUIRE (FILE='image_size.dat', exist=ex_IM)

         if (ex_im ) then
         open(1,file=    'image_size.dat')
         DO I=1,10000
            read(1,*,end=22) timeim(I),rIM(I)
         enddo
         
 22      NIM=I-1
      continue

      endif
      endif
      rcore=0d0

      if (ex_im ) then
         else
            return
         endif
         
      
      n1=1
      
         x1(1) = time
         
         call intpl(timeIM, rIM, nIM, x1,y1,n1)
         r1=y1(1)
         
      pi=dacos(-1d0)
        vliq=2*pi*rliq**3/3

        rcore=0d0
        

          
          if ( r1.gt.rliq) then


             rm=r1
             xmin=rliq
             xmax=rm
             erabs=0d0
             errel=0d0
             itmax=100
             call dzbren(fvol, erabs,errel,xmin,xmax,ITMAX)
             
             ra= xmax
             rcore = (ra**3-rliq**3)**(1/3d0)
c             print *, 'ra, f', ra, fvol(ra), rcore
             if (rcore.le.0)rcore=0d0

             
c            call getcore(vliq,rkappe,rdrop,rsphere,rcore,thick)
c        print*, r1, rliq, rcore
            
            
         endif
c
c         print*,'aa'
         
c         stop
         
         
c         write(6,'(4E15.6)') time1, r1, rliq,rcore

c     recalculate rkappe from rcore and rsphere


         return
         end
      
      
c     v: volume of the droplet (Liquid +solid, exclusice core)
c     rkappe: radius of image
      

      

      function fvol(ra)
      implicit real*8 (a-h,o-z)
      common /fvoll/ rm, rliq
      fvol=ra**3-( (rm-ra)/(dsqrt(2d0)-1))**3 -rliq**3
c      print*, 'f ' , ra**3-( (rm-ra)/(dsqrt(2d0)-1))**3 , rliq**3,rliq
      
      return
      end
      

c     importtant
      subroutine        caldl(ta,aw,ml,x,J,dl)
      implicit real*8 (a-h,o-z)
      parameter (np=50, nsmm=100)
      integer izc(np)
      real*8 mm(np),mv(np),x(*)
      
      common /isco2/iseqco2,is,nss,isupdate,iseqNH3,iseqaa
      
      common /M/ MM,mv,izc
      real*8 ml(np)
      real*8 xn(NSMM,NP)
      common /xn/ xn

      common /xvion/ xvion,dl_ion,dl_ions(NSMM)
      common /output/imode_output,idiff,imode_pH,imode_eq,imode_NH4NO3
      
      
c     aw0=awshell(J)

c         ienh=0
c     enhanced diff only when solid is in the center
c         if (j.eq.16 .and. xnsolid.gt.0d0 .and. iscenter .eq.1) ienh=1
c         if (j.eq.17 .and. xnsolid.gt.0d0 .and. iscenter .eq.1) ienh=1

         
      if (idiff.eq.0) then
         xvv= mv(j)
      if( izc(j).ne.0) xvv=xvion

       call cal_dlaw_eta(Ta,aw,dl,xvv)
      endif
      
c       call cal_dlaw(Ta,aw,dl)  !diffusin coefficient of ions
c        call cal_dlaw_walker_mod(Ta,aw,dl)

      if (idiff.eq.3)   call cal_dlaw_suc(Ta,aw,dl) 
      if (idiff.eq.4)   call cal_dlaw_citric(Ta,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(Ta,aw,dl)
      if (idiff.eq.6)   call cal_dlaw_walker_h2o(Ta,aw,dl)
       

      if( izc(j).eq.0 ) then
      if (idiff.eq.5)   call cal_dlaw_walker(Ta,aw,dl)
       endif

       if ( idiff.eq.10 )    then



c     Dl of ions
                  call DH2O_NaCO3_cl(Ta,aw,ml,dl1,dlion)
c     Dl of neutral species
                  call DH2O_Nacl_suc(Ta,ML,aw,dl)
c                  print*,'d10',dl ,dlion
                  
                  
                  if( (izc(j)).ne.0) then
                 dl=dlion
                endif
         endif

       
       if ( idiff.eq.12 )    then

            xm0=ml(16)

                  call DH2O_NaCO3_cl(Ta,aw,ml,dl,dlion)
      if( (izc(j)).ne. 0)         dl=dlion


         endif

      if (idiff.eq.13) then
c     take the mean value
            xmv=mv(j)


             if (izc(j).ne.0) xmv=xvion
              call cal_dl_an_suc(Ta,ML,aw,xmv,dl)           
              ratio = ml(18)/(ml(2)+ml(18))
              if (dl_ion.gt.dl) dl_ion=dl
              if (izc(j).ne.0 .and. imode_MA.ne.1 ) dl_ions(I)=dl !ions
              if (imode_MA.eq.1.and. j.eq.29) dl_ions(I)=dl  !Maleic acid
           endif
       
c
c     set uplimit for EDB 
            if (imode_output.eq.2 .and. imode_NH4NO3+imode_MA.ge.1) then
               if (is.le.0) is=1
               if (is.gt.nss) is=NSs
               
               dl_upper = (x(Is+1))**2/10 !  10s diffusin length
                 if (dl.ge.dl_upper) dl= dl_upper
              endif
              
              return
              end
              
c
c     calculate the diffusion coeffcieint of H2O as a function of aw obatined from from present EDB study used in present paper
c
c     ----------------------------------------------------------------------------------------
   
      subroutine cal_dlaw_walker_mod(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)

      aw1=1

      t0=298.15
      call cal_dlaw_citric(T0,aw1,dl298)
      t293=293.15
      call cal_dlaw_citric(T293,aw1,dl293)
c      print*, dl0
      
      Dl0= dl293 * 2.44D-5/dl298  ! scale d0 to 2.44E-5 cm2/s at 298d0 and aw=1
           
            fac(1)=7E-9*5
            fac(2)=8E-9*5   !10E-9
           fac(3)= 8D-9*5 !11D-9

           fac(4)= 9E-9*5

       fac(5)= 1.8E-8*5
       fac(6)= 6E-8*5
       fac(7)= dl0/4
       fac(8)= dl0*.8
       fac(9)= dl0

      rh(1)= 0
      rh(2)=.5
       rh(3)=.65
       rh(4)=.75

       rh(5)=0.85
       rh(6)=0.90
       rh(7)=0.96
       rh(8)=0.99
       
       rh(9)=1

      fac=dlog(fac)
       N1=1
       N5=9
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))

       t0=293.15
       call  cal_dlaw_citric(t,aw0,dlt)
       call  cal_dlaw_citric(t0,aw0,dlt0)
       dl= dl* dlt/dlt0
       
      

       return
      end
