
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

      integer                       :: MPI_PART_XFT
      !type(mpi_part_t), allocatable :: part_portal(:)
      !type(mpi_part_t), allocatable :: sbuff(:),rbuff(:) ! don't need to be targets?

      type(mpi_part_t), dimension(:), allocatable :: sbuf,rbuf


      integer :: iblk(2),idtp(2) 
      INTEGER(KIND=MPI_ADDRESS_KIND) :: idis(2)

      TYPE(mpi_part_t) :: pdum(2)

      integer :: k,r,s,i,itag,m
      integer :: icomw,id,ier,nproc
      integer, dimension(:), allocatable :: scounts,rcounts
      integer :: mpistat(MPI_STATUS_SIZE)
      !integer, dimension(:), allocatable :: ir

      integer :: my_nsend,my_nrecv
      integer, dimension(:), allocatable :: my_ps,my_pr

      type(part_t), pointer :: pPtr
      type(part_t) :: item

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
      if (nproc.ne.4) STOP

      ! send and receive buffer, oversized for now
      ! eventually fix size and do while particles remain to send or receive
      ! i.e. "particle packets"
      !allocate(sbuff(MAX_NUM),rbuff(MAX_NUM))

      scounts = 0
      rcounts = 0
 
      if (id.eq.0) then
         scounts(1) = 3 ! send 3 particles to rank 1
      endif

      ! START
      
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

      if (id.eq.0) then

         do i=1,3
            item%pid = i
            call insert_item(item)
         enddo

      endif

      !call print_list()

      do k=1,my_nsend
         m = my_ps(k)

         itag = 1000 + m

#ifdef _DEBUG_
         PRINT*,'rank',id,'sending',scounts(m),'particles to rank',m
#endif

         if (allocated(sbuf)) deallocate(sbuf)
         allocate(sbuf(scounts(m)))

         ! copy from link-list to send buffer
         ! loop though linked-list, copy items, then delete
         ! eventually do it item-by-item

         i = 0
         if (associated(head)) then
            pPtr => head
            do while (associated(pPtr))
               i = i + 1
               sbuf(i)%pid = pPtr%pid
               pPtr => pPtr%next
            enddo
         endif
         if (i.ne.scounts(m)) STOP

         call MPI_SEND(sbuf(1),scounts(m),MPI_PART_XFT, &
                       m,itag,icomw,ier)

         do i=1,scounts(m)
            PRINT*,sbuf(i)%pid
         enddo 

      enddo

      do k=1,my_nrecv
         m = my_pr(k)

         itag = 1000 + id

#ifdef _DEBUG_
         PRINT*,'rank',id,'receiving',rcounts(m),'particles from rank',m 
#endif

         if (allocated(rbuf)) deallocate(rbuf)
         allocate(rbuf(rcounts(m))) 

         call MPI_RECV(rbuf(1),rcounts(m),MPI_PART_XFT, &
                       m,itag,icomw,mpistat,ier)

         do i=1,rcounts(m)
            PRINT*,rbuf(i)%pid
         enddo 

      enddo

      ! here I know I am only sending from 0 to 1
      ! but eventually I will need to break up the
      ! linked-list and copy to the send buffer
      ! in packets, for different ranks 

      call MPI_Finalize(ier)
    
      end program


























!      ict = sum(scounts) + sum(rcounts)
!
!      ! do while particles remain to send or receive
!      do while (ict.gt.0)
! 
!         ! loop over peers to send, m
!         do n=1,n_peers 
!   
!            m = my_peers(n)
!
!            ! if (communication of packet to peer m was completed
             !     AND I still have particles to send to m (scount(m).gt.0?) ) then 

!            ! copy data from linked-list to send buffer for peer m
!
!            ! TESTING
!            sbuff(:)%pid = 0 
!            do i=1,bsize !scounts(m)
!               sbuff(i)%pid = i
!            enddo
!
!            itag = 1000 + m ! constant + id of receiver
!
!            !if (id.eq.0) PRINT*,'sending from ' 
!
!            call MPI_SEND(sbuff(1),bsize,MPI_PART_XFT, &
!                          m,itag,icomw,ier)
!
!            ! reduce send count
!            scounts(m) = max(0,scounts(m)-bsize)
! 
!            ! immediately clear send buffer for testing
!            sbuff(:)%pid = 0
!   
!   !         call MPI_ISEND(part_portal(1),scount(m),MPI_PART_XFT, &
!   !                        m,itag,icomw,ir(m),ier)
!   
!         enddo
!   
!      ! loop over peers from which to receive, n
!
!      do n=1,n_peers
!
!         m = my_peers(n)
!
!         itag = 1000 + id ! constant + id of receiver (me)
!
!
!         call MPI_RECV(rbuff(1),bsize,MPI_PART_XFT, &
!                       m,itag,icomw,mpistat,ier)
!
!         
!!         if (rcounts(m).gt.0) then
!!            PRINT*,'rank',id,'received the following',rcounts(m),'items from rank',m,':'
!!            !do i=1,rcounts(m)
!!            !   PRINT*,rbuff(i)%pid
!!            !enddo
!!         endif
!
!         rcounts(m) = max(0,rcounts(m)-bsize)
!
!      enddo
!
!         ! recompute total count
!         ict = sum(scounts) + sum(rcounts)
!
!      enddo     
