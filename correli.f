      parameter(nmax=4096)
      DIMENSION sismo1(4096),sismo2(4096),res(4096),x(4096),
     $          vec1(4096),vec2(4096)
      complex cor(8192)
      integer sizear
      character*100 arch1,arch2,outpfile,outpfile1,arch1n,arch2n
      write(*,*)'File no 1'
      read(*,*)arch1
      write(*,*)'File no 2'
      read(*,*)arch2
      sizear=lnblnk(arch1)
      arch1n=arch1(1:sizear-3)//'dat'
      write(*,*)arch1n
      write(*,*)sizear
      sizear=lnblnk(arch2)
      arch2n=arch2(1:sizear-3)//'dat'
      write(*,*)arch2n
      write(*,*)sizear
      open (16,file=arch1n)
      open (17,file=arch2n)

c     arch1='uno.1.sac'
c     arch2='cuatro.1.sac'
      ventana=0.3
      call RSAC1(arch1,sismo1,npts1,beg,delta,4096,nerr)
      if (nerr.ne.0) then
          write(*,*)'Some problem reading file ',arch1
          write(*,*)nerr
          stop
      endif
      call RSAC1(arch2,sismo2,npts2,beg,delta,4096,nerr)
      if (nerr.ne.0) then
          write(*,*)'Some problem reading file ',arch2
          stop
      endif
c     Renaming files
      outpfile = 'CORR.'//arch1(45:46)//'-'//arch2(45:46)//'.'
     $           //arch1(5:15)//'.dat'
      write(*,*)'Is the following name ok?, (if yes just press enter)'
      write(*,*)outpfile
      read(*,'(a)')outpfile1
      if(outpfile1(1:1).eq.'')outpfile1=outpfile
      open(15,file = outpfile1)
      limite1=1
      limite2=int(ventana/delta)
      nptsmax=npts1
      if(npts1.ge.npts2)nptsmax=npts2
      ancdelta=0.0
      aintervalo=0.0
      do j=1,4000
         if (limite1.gt.nptsmax.or.limite2.gt.nptsmax)goto 1999
c        Define time window
         aintervalo=aintervalo+ventana
         limite2=limite2+1
         k=0
         do i=1,nptsmax
            vec1(i)=0.0
            vec2(i)=0.0
         enddo
         aint1=0.0
         aint2=0.0
         do i=limite1,limite2
            k=k+1
            vec1(k)=sismo1(i)
            aint1=aint1+abs(sismo1(i))**2
            vec2(k)=sismo2(i)
            aint2=aint2+abs(sismo2(i))**2
         enddo
c        limite1=limite1+1
          
         if (npts1.gt.4096)then
            write(*,*)'exiting because number of points'
            goto 1999
         endif
         if (npts1.le.4096.or.npts2.le.1024)npts=4096
         if (npts1.le.2048.or.npts2.le.1024)npts=2048
         if (npts1.le.1024.or.npts2.le.1024)npts=1024
         if (npts1.le.512.or.npts2.le.512)npts=512
         if (npts1.le.256.or.npts2.le.256)npts=256
         if (npts1.le.128.or.npts2.le.128)npts=128
         if (npts1.le.64.or.npts2.le.64)npts=64
         do i=1,npts
            res(i)=0.0
         enddo
          
         call correl(vec2,vec1,npts,cor)

         k=(npts/4)+1
         do i=1,npts/4
            res(k)=real(cor(i))
            k=k+1
         enddo
         k=(npts/4)+1
c        write(*,*)k
         do i=1,npts/4
            res(i)=real(cor(k))
            k=k+1
         enddo
         corrim=-(delta*npts)/2.0
         delta1=delta*2.0
*        write(*,*)limite1,limite2,corrim,delta1,npts/4,nptsmax,
*    $             npts1,npts2
c        call wsac1('salida.sac',res,npts/2,corrim,delta1,nerr)
c        call wsac1('sismo1.sac',vec1,npts/2,corrim,delta,nerr)
c        call wsac1('sismo2.sac',vec2,npts/2,corrim,delta,nerr)
c        Calculando la integral para normalizacion
         fac=(aint1*aint2)**0.5
c        Buscando el maximo
         amax=0.0000
         do i=1,nptsmax
              if(res(i).gt.amax)amax=res(i)
         enddo
         amaxc=amax/fac
c        Performing noise correction
         ll=0
         u1=0.0
         u2=0.0
         an1=0.
         an2=0.
         do i=limite1,limite2
            ll=ll+1
            u1=u1+vec1(i)**2
            u2=u1+vec2(i)**2
         enddo
            u1=u1/float(ll)
            u2=u2/float(ll)
         if (j.eq.1) then
            ll=0
            do i=limite1,limite2
               ll=ll+1
               an1=an1+vec1(i)**2
               an2=an1+vec2(i)**2
            enddo
            an1=an1/float(ll)
            an2=an2/float(ll)
         endif
         anoisec=(sqrt(1-(an1/u1))*sqrt(1-(an2/u2)))
         amaxc=amaxc/anoisec
         if (amaxc.gt.1.0)amaxc=1.0
c        Considerando lp 5 hz (promedio 3hz)
c         alpha=3897.0  Valores originales
c         beta=2250.0   Valores originales
         alpha=8313.0
         beta=4800.0
         sigma=2*(1-amaxc)/355.3
         anum=7*((2/alpha**6)+(3/beta**6))
         aden=(6/alpha**8)+(7/beta**8)
         dist=sqrt(sigma*anum/aden)
         write(15,15)float(limite2-limite1)*delta,amaxc,dist
         ancdelta=ancdelta+delta
         write(16,*)ancdelta,sismo1(j)
         write(17,*)ancdelta,sismo2(j)
      enddo
1999  continue
14    format(i4,6f8.3)
15    format(3f8.3)
      end
