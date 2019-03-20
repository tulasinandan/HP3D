!******************************************************************************
!             3D particle code: random number generator
!                             Marc Swisdak, 2007
!
!                        LATEST CHANGE: October 1, 2007
!
!******************************************************************************

! These routines are shamelessly copied from _Numerical Recipes in 
! Fortran 90_, Chapter B7.  Earlier versions of p3d used the F90 
! intrinsic calls to random_seed() and random_number(), but not
! surprisingly those proved to be platform-dependent.  Eventually
! p3d was ported to a machine (franklin) on which those intrinsic
! functions were implemented badly.  Hence the need for including
! our own.

! This file includes the support module ran_state as well as the
! generator ran1.

MODULE ran_state
	USE nrtype
	IMPLICIT NONE
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
	INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
	INTEGER(K4B), SAVE :: lenran=0, seq=0
	INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
	INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
	INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
		nran,mran,ranv
	REAL(SP), SAVE :: amm
	INTERFACE ran_hash
		MODULE PROCEDURE ran_hash_s, ran_hash_v
	END INTERFACE
CONTAINS
!BL
	SUBROUTINE ran_init(length)
	USE nrtype; USE nrutil, ONLY : arth,nrerror,reallocate
	IMPLICIT NONE
	INTEGER(K4B), INTENT(IN) :: length
	INTEGER(K4B) :: new,j,hgt
	if (length < lenran) RETURN
	hgt=hg
	if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
	if (hgng >= 0)        call nrerror('ran_init: arith assump 2 fails')
	if (hgt+1 /= hgng)    call nrerror('ran_init: arith assump 3 fails')
	if (not(hg) >= 0)     call nrerror('ran_init: arith assump 4 fails')
	if (not(hgng) < 0)    call nrerror('ran_init: arith assump 5 fails')
	if (hg+hgng >= 0)     call nrerror('ran_init: arith assump 6 fails')
	if (not(-1_k4b) < 0)  call nrerror('ran_init: arith assump 7 fails')
	if (not(0_k4b) >= 0)  call nrerror('ran_init: arith assump 8 fails')
	if (not(1_k4b) >= 0)  call nrerror('ran_init: arith assump 9 fails')
	if (lenran > 0) then
		ranseeds=>reallocate(ranseeds,length,5)
		ranv=>reallocate(ranv,length-1)
		new=lenran+1
	else
		allocate(ranseeds(length,5))
		allocate(ranv(length-1))
		new=1
		amm=nearest(1.0_sp,-1.0_sp)/hgng
		if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
			call nrerror('ran_init: arth assump 10 fails')
	end if
	ranseeds(new:,1)=seq
	ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
	do j=1,4
		call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
	end do
	where (ranseeds(new:,1:3) < 0) &
		ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
	where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
	if (new == 1) then
		iran0=ranseeds(1,1)
		jran0=ranseeds(1,2)
		kran0=ranseeds(1,3)
		mran0=ranseeds(1,4)
		nran0=ranseeds(1,5)
		rans=nran0
	end if
	if (length > 1) then
		iran => ranseeds(2:,1)
		jran => ranseeds(2:,2)
		kran => ranseeds(2:,3)
		mran => ranseeds(2:,4)
		nran => ranseeds(2:,5)
		ranv = nran
	end if
	lenran=length
	END SUBROUTINE ran_init
!BL
	SUBROUTINE ran_deallocate
	if (lenran > 0) then
		deallocate(ranseeds,ranv)
		nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
		lenran = 0
	end if
	END SUBROUTINE ran_deallocate
!BL
	SUBROUTINE ran_seed(sequence,size,put,get)
	IMPLICIT NONE
	INTEGER, OPTIONAL, INTENT(IN) :: sequence
	INTEGER, OPTIONAL, INTENT(OUT) :: size
	INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
	INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
	if (present(size)) then
		size=5*lenran
	else if (present(put)) then
		if (lenran == 0) RETURN
		ranseeds=reshape(put,shape(ranseeds))
		where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
		where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
		iran0=ranseeds(1,1)
		jran0=ranseeds(1,2)
		kran0=ranseeds(1,3)
		mran0=ranseeds(1,4)
		nran0=ranseeds(1,5)
	else if (present(get)) then
		if (lenran == 0) RETURN
		ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
		get=reshape(ranseeds,shape(get))
	else if (present(sequence)) then
		call ran_deallocate
		seq=sequence
	end if
	END SUBROUTINE ran_seed
!BL
	SUBROUTINE ran_hash_s(il,ir)
	IMPLICIT NONE
	INTEGER(K4B), INTENT(INOUT) :: il,ir
	INTEGER(K4B) :: is,j
	do j=1,4
		is=ir
		ir=ieor(ir,ishft(ir,5))+1422217823
		ir=ieor(ir,ishft(ir,-16))+1842055030
		ir=ieor(ir,ishft(ir,9))+80567781
		ir=ieor(il,ir)
		il=is
	end do
	END SUBROUTINE ran_hash_s
!BL
	SUBROUTINE ran_hash_v(il,ir)
	IMPLICIT NONE
	INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
	INTEGER(K4B), DIMENSION(size(il)) :: is
	INTEGER(K4B) :: j
	do j=1,4
		is=ir
		ir=ieor(ir,ishft(ir,5))+1422217823
		ir=ieor(ir,ishft(ir,-16))+1842055030
		ir=ieor(ir,ishft(ir,9))+80567781
		ir=ieor(il,ir)
		il=is
	end do
	END SUBROUTINE ran_hash_v
END MODULE ran_state


	SUBROUTINE ran1_s(harvest)
	USE nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran0,jran0,kran0,nran0,mran0,rans
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	if (lenran < 1) call ran_init(1)
	rans=iran0-kran0
	if (rans < 0) rans=rans+2147483579_k4b
	iran0=jran0
	jran0=kran0
	kran0=rans
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	if (nran0 == 1) nran0=270369_k4b
	mran0=ieor(mran0,ishft(mran0,5))
	mran0=ieor(mran0,ishft(mran0,-13))
	mran0=ieor(mran0,ishft(mran0,6))
	rans=ieor(nran0,rans)+mran0
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran1_s

	SUBROUTINE ran1_v(harvest)
	USE nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran,jran,kran,nran,mran,ranv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	ranv(1:n)=iran(1:n)-kran(1:n)
	where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
	iran(1:n)=jran(1:n)
	jran(1:n)=kran(1:n)
	kran(1:n)=ranv(1:n)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	where (nran(1:n) == 1) nran(1:n)=270369_k4b
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
	ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran1_v
