
      ! linked-list format
!      type, public :: part_t
!         integer :: pid = 0       ! (global) particle index
!         integer :: i = 0         ! local cell index that contains the particle
!         real*8 :: xyz(3) = 0.0d0 ! particle position
!         real*8 :: uvw(3) = 0.0d0 ! particle velocity
!         real*8 :: tleft  = 0.0d0 ! remaining time until particle has moved dt
!         type(part_t), pointer :: next=>NULL()
!         !type(part_t), pointer :: prev=>NULL()
!      end type

      program main

      use llist

      use mpi

      implicit none

      !integer, parameter :: bsize = 2 ! fixed buffer size

      type :: mpi_part_t
         integer :: pid = 0        ! particle index
         integer :: i = 0          ! local cell index that contains the particle
         real*8  :: xyz(3) = 0.0d0 ! particle position
         real*8  :: uvw(3) = 0.0d0 ! particle velocity
         real*8  :: tleft  = 0.0d0 ! ! remaining time until particle has moved dt
      end type

      integer :: MPI_PART_XFT

      type(mpi_part_t), dimension(:), allocatable :: sbuf,rbuf

      integer :: iblk(2),idtp(2) 
      integer(kind=MPI_ADDRESS_KIND) :: idis(2)

      type(mpi_part_t) :: pdum(2)

      integer :: k,r,s,i,itag,m
      integer :: icomw,id,ier,nproc
      integer, dimension(:), allocatable :: scounts,rcounts
      integer :: mpistat(MPI_STATUS_SIZE)
      !integer, dimension(:), allocatable :: ir

      integer :: my_nsend,my_nrecv
      integer, dimension(:), allocatable :: my_ps,my_pr
      integer, dimension(:), allocatable :: sdis,rdis

      type(part_t), pointer :: pPtr
      type(part_t) :: item

      !integer :: ib,ie

      call MPI_Init(ier)

      icomw = MPI_COMM_WORLD

      call MPI_Comm_rank(icomw,id,ier)
      call MPI_Comm_size(icomw,nproc,ier)

      allocate(scounts(0:nproc-1),rcounts(0:nproc-1))
      !allocate(ir(0:2*nproc-1))

      iblk = (/2, 7/) ! 2 integers, 7 (2*3+1) reals
     
      CALL  MPI_GET_ADDRESS(pdum(1)%pid    ,idis(1),ier) ! first integer
      CALL  MPI_GET_ADDRESS(pdum(1)%xyz(1) ,idis(2),ier) ! first real
      idis(2) = idis(2) - idis(1)  ! offset until reals
      idis(1) = 0

      idtp(1) = MPI_INTEGER
      idtp(2) = MPI_DOUBLE_PRECISION

      CALL MPI_TYPE_CREATE_STRUCT(2, iblk, idis, idtp, MPI_PART_XFT, ier)
      CALL MPI_TYPE_COMMIT(MPI_PART_XFT,ier)

      ! four processor example for now
      !if (nproc.ne.4) STOP

      ! send and receive buffer, oversized for now
      ! eventually fix size and do while particles remain to send or receive
      ! i.e. "particle packets"

      scounts = 0
      rcounts = 0
 
      if (id.eq.0) then
         scounts(1) = 3 ! send 3 particles to rank 1
         scounts(3) = 9
         scounts(6) = 2
      endif

      ! FILL DUMMY LINKED-LIST(S)
      do i=1,sum(scounts)
         item%pid = i
         call insert_item(item)
      enddo
      !call print_list()

      allocate(sbuf(sum(scounts)))


 
      ! improve by only sending/receiving
      ! between partition shared ranks
      call MPI_ALLTOALL(scounts,1,MPI_INTEGER, &
                        rcounts,1,MPI_INTEGER, &
                        icomw,ier)

      ! resize these every time of leave sized to
      ! max. number of possible peers and allow for a 
      ! send of size zero (which is likely not that common)?
      ! in this way, we can use in_nsend/in_pr

      my_nsend = 0
      my_nrecv = 0

      do m=0,nproc-1
         if (scounts(m).gt.0) my_nsend = my_nsend + 1
         if (rcounts(m).gt.0) my_nrecv = my_nrecv + 1
      enddo

      if (allocated(my_ps)) deallocate(my_ps)
      if (allocated(my_pr)) deallocate(my_pr)

      allocate(my_ps(my_nsend)) ! index from 1
      allocate(my_pr(my_nrecv)) ! index from 1

      s = 1 ! index from 1
      r = 1 ! index from 1

      do m=0,nproc-1
         if (scounts(m).gt.0) then
            my_ps(s) = m
            s = s + 1
         endif
          if (rcounts(m).gt.0) then
            my_pr(r) = m
            r = r + 1
         endif
      enddo

!      allocate(sdis(my_nsend))
!      allocate(rdis(my_nrecv))

      ! create displacement arrays (if there is anything to send/recv)
      if (my_nsend>0) then
         allocate(sdis(my_nsend))
         sdis(1) = 0 
         do k=2,my_nsend
            sdis(k) = sdis(k-1) + scounts(my_ps(k-1))
         enddo
      endif

      if (my_nrecv>0) then
         allocate(rdis(my_nrecv))
         rdis(1) = 0
         do k=2,my_nrecv
            rdis(k) = rdis(k-1) + rcounts(my_pr(k-1))
         enddo
      endif

!      do k=1,my_nsend
!         m = my_ps(k)
!         ib = sdis(k)+1
!         ie = ib + scounts(m) - 1
!      enddo

      ! start at the head of the linked-list
      if (associated(head)) pPtr => head

      do k=1,my_nsend
         m = my_ps(k)

         itag = 1000 + m

#ifdef _DEBUG_
         PRINT*,'rank',id,'sending',scounts(m),'particles to rank',m
#endif

         ! copy from link-list to send buffer
         ! loop though linked-list, copy items, then delete
         ! eventually do it item-by-item

         ! traverse the linked-list
         do i=1,scounts(m)
            sbuf(sdis(k)+i)%pid = pPtr%pid
            pPtr => pPtr%next
         enddo

!         if (associated(head)) then
!            pPtr => head
!            do while (associated(pPtr))
!               i = i + 1
!               sbuf(i)%pid = pPtr%pid
!               pPtr => pPtr%next
!            enddo
!         endif

         call MPI_SEND(sbuf(sdis(k)+1),scounts(m),MPI_PART_XFT, &
                       m,itag,icomw,ier)

         !do i=1,scounts(m)
         !   PRINT*,sbuf(i)%pid
         !enddo 

      enddo

      !allocate(rbuf(maxval(rcounts))) ! re-use

      do k=1,my_nrecv
         m = my_pr(k)

         itag = 1000 + id

#ifdef _DEBUG_
         PRINT*,'rank',id,'receiving',rcounts(m),'particles from rank',m 
#endif

         if (allocated(rbuf)) deallocate(rbuf)
         allocate(rbuf(rcounts(m)))

         call MPI_RECV(rbuf(rdis(k)+1),rcounts(m),MPI_PART_XFT, &
                       m,itag,icomw,mpistat,ier)

         do i=1,rcounts(m)
            PRINT*,id,'|',rbuf(i)%pid
         enddo 

      enddo

      ! here I know I am only sending from 0 to 1
      ! but eventually I will need to break up the
      ! linked-list and copy to the send buffer
      ! in packets, for different ranks 

      call MPI_Finalize(ier)
    
      end program

