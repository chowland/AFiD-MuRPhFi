    forclo = 1.d0
    if(infig.gt.0) then
        usaldto = 1./aldto
        do n=1,npunx
            ic=indgeo(3,n,1)
            jc=indgeo(3,n,2)
            kc=indgeo(3,n,3)
            forclo(kc,jc,ic)=0.d0
            !           ie=indgeoe(3,n,1)
            !           je=indgeoe(3,n,2)
            ke=indgeoe(3,n,3)
            q3e=((al*dt+aldto)*vx(ke,jc,ic)-al*dt*q3bo(n))*usaldto
            rhs(kc,jc,ic) = -vx(kc,jc,ic) + q3e*distb(3,n)
            q3bo(n)= vx(ke,jc,ic)
        end do
    end if
    
    call SolveImpEqnUpdate_X_ibm
