       forclo=1.d0 
       if(infig.gt.0) then
         usaldto = 1./aldto
         do n=1,npunte
           ic=indgeot(n,1)
           jc=indgeot(n,2)
           kc=indgeot(n,3)
           forclo(kc,jc,ic)=0.d0
!           ie=indgeoet(n,1)
!           je=indgeoet(n,2)
           ke=indgeoet(n,3)
           dense=((al*dt+aldto)*temp(ke,jc,ic)-al*dt*densb(n))*usaldto
           rhs(kc,jc,ic) = -temp(kc,jc,ic) + dense*distbt(n)  &
                         +(1.-distbt(n))*temb(n)
           densb(n)= temp(ke,jc,ic)

           if(ifnoslipy.eq.1) then     ! SL
               if (jc==1) then
                    rhs(kc,jc,ic)=temp(kc,jc+1,ic)-temp(kc,jc,ic)
               elseif (jc==nym) then
                    rhs(kc,jc,ic)=temp(kc,jc-1,ic)-temp(kc,jc,ic)
               endif
           endif
           
           if(ifnoslipz.eq.1) then     ! SL
               if (ic==1) then
                    rhs(kc,jc,ic)=temp(kc,jc,ic+1)-temp(kc,jc,ic)
               elseif (ic==nzm) then
                    rhs(kc,jc,ic)=temp(kc,jc,ic-1)-temp(kc,jc,ic)
               endif
           endif
           
         end do
        end if
       
         
         call SolveImpEqnUpdate_Temp_ibm
