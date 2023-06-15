!RO    Extra instructions in the case of Immersed boundary for invtr1
       forclo=1.d0

       if(infig.gt.0) then
         usaldto = 1./aldto
         do n=1,npunz
           ic=indgeo(1,n,1)
           jc=indgeo(1,n,2)
           kc=indgeo(1,n,3)
           forclo(kc,jc,ic)=0.d0
!           ie=indgeoe(1,n,1)
!           je=indgeoe(1,n,2)
           ke=indgeoe(1,n,3)
           q1e=((al*dt+aldto)*vz(ke,jc,ic)-al*dt*q1bo(n))*usaldto
           rhs(kc,jc,ic) = -vz(kc,jc,ic) + q1e*distb(1,n)
           q1bo(n)= vz(ke,jc,ic)
         end do
        endif
        
        
        call SolveImpEqnUpdate_YZ_ibm(vz,rhs,forclo)
        
