!RO    Extra instructions in the case of Immersed boundary for invtr2

       forclo=1.d0
       if(infig.gt.0) then
         usaldto = 1./aldto
         do n=1,npuny
             ic=indgeo(2,n,1)
             jc=indgeo(2,n,2)
             kc=indgeo(2,n,3)
             forclo(kc,jc,ic)=0.d0
!             ie=indgeoe(2,n,1)
!             je=indgeoe(2,n,2)
             ke=indgeoe(2,n,3)
             q2e=((al*dt+aldto)*vy(ke,jc,ic)-al*dt*q2bo(n))*usaldto
             rhs(kc,jc,ic) = -vy(kc,jc,ic) + q2e*distb(2,n)
             q2bo(n)= vy(ke,jc,ic)
         end do
        endif
         
        call SolveImpEqnUpdate_YZ_ibm(vy,rhs,forclo)
       
