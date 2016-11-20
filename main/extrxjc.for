      module MOD_EXTRXJC
      contains      
      subroutine extrxjc (XJCREA,XJC,EDDREA,EDDI,nd,nf,K)
c***     now extract XJC and EDDI for the frequency K
      IMPLICIT REAL*8(A-H,O-Z)

      dimension xjcrea(nd,nf),xjc(nd),eddrea(3,nd,nf),eddi(3,nd)
	! xjc(:) = xjcrea(:,K)
	! eddi(1:3,:) = eddrea(1:3,:,K)
	do L=1,nd
	   xjc(L)=xjcrea(L,K)
	   eddi(1,L)=eddrea(1,L,K)
	   eddi(2,L)=eddrea(2,L,K)
	   eddi(3,L)=eddrea(3,L,K)
      enddo

	return
	end subroutine
      end module