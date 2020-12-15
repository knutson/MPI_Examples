
      ! eventually fix buffer size and do while particles remain to send or receive i.e. "particle packets"

      program main

      use llist

      use mpi

      implicit none

      type :: mpi_part_t
         ! pid and id should not need to be exchanged       
         integer :: i       = -1          ! local index of cell that contains the particle
         real*8 :: xyz(3)   = huge(1.0d0) ! particle position
         real*8 :: uvw(3)   = huge(1.0d0) ! current particle velocity
         real*8 :: uvw_n(3) = huge(1.0d0) ! previous particle velocity
         real*8 :: tleft    = huge(1.0d0) ! remaining time until particle has moved pdt (particle dt) 
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

      integer :: pid

      call MPI_Init(ier)

      icomw = MPI_COMM_WORLD

      call MPI_Comm_rank(icomw,id,ier)
      call MPI_Comm_size(icomw,nproc,ier)

      allocate(scounts(0:nproc-1),rcounts(0:nproc-1)) ! over-sized for now
      !allocate(ir(0:2*nproc-1))

      iblk = (/1, 10/) ! 1 integers, 10 (3*3+1) reals
     
      CALL  MPI_GET_ADDRESS(pdum(1)%i      ,idis(1),ier) ! first integer
      CALL  MPI_GET_ADDRESS(pdum(1)%xyz(1) ,idis(2),ier) ! first real
      idis(2) = idis(2) - idis(1) ! offset until reals
      idis(1) = 0                 ! offset until integers (zero)

      idtp(1) = MPI_INTEGER
      idtp(2) = MPI_DOUBLE_PRECISION

      CALL MPI_TYPE_CREATE_STRUCT(2, iblk, idis, idtp, MPI_PART_XFT, ier)
      CALL MPI_TYPE_COMMIT(MPI_PART_XFT,ier)
 
      ! initialize the hat (this is where PIDs are set)
      call init_hat()

!      ! initialize linked-list
!      if (id.eq.0) then
!         do i=1,14
!            pid = pull_from_hat()
!            pPtr => particles(pid)%p
!            call insert_link(pPtr)
!#ifdef _DEBUG_
!            ! at this point, i should equal pid
!            if (i.ne.pid) then
!               PRINT*,'Error: i = ',i,', PID = ',pid
!               STOP
!            endif
!#endif
!            !PRINT*,'particle count on id',id,'=',npart
!         enddo
!      endif

      ! fill the linked-list (copy from static array or buffer)
      if (id.eq.0) then
         do i=1,1 ! send one to self
            pid = pull_from_hat()
            pPtr => particles(pid)%p
            call insert_link(pPtr)
            !pPtr => particles(i)%p
            pPtr%id = 0
            pPtr%tleft = dble(i)
         enddo
         do i=2,4
            pid = pull_from_hat()
            pPtr => particles(pid)%p
            call insert_link(pPtr)
            !pPtr => particles(i)%p
            pPtr%id = 1
            pPtr%tleft = dble(i)
         enddo
         do i=5,12
            pid = pull_from_hat()
            pPtr => particles(pid)%p
            call insert_link(pPtr)
            !pPtr => particles(i)%p
            pPtr%id = 3
            pPtr%tleft = dble(i)
         enddo
         do i=13,14
            pid = pull_from_hat()
            pPtr => particles(pid)%p
            call insert_link(pPtr)
            !pPtr => particles(i)%p
            pPtr%id = 6
            pPtr%tleft = dble(i)
         enddo
      endif

      ! initialize and fill linked-lists
!      if (id.eq.1) then
!         do i=1,5
!            pid = pull_from_hat()
!            pPtr => particles(pid)%p
!            call insert_link(pPtr)
!            !if (i.ne.pid) STOP ! relax this constraint
!            !pPtr => particles(i)%p
!
!            pPtr%id = 0
!            pPtr%tleft = dble(i)
!         enddo
!      endif

      !call print_list()
      !PRINT*,'size of linked-list owned by id',id,'=',get_count()

      scounts = 0
      rcounts = 0

      ! determine send counts by traversing linked-list and counting
      if (associated(head)) then
         pPtr => head
         do while (associated(pPtr))
            !!!if (pPtr%id.ne.id) scounts(pPtr%id) = scounts(pPtr%id) + 1
            scounts(pPtr%id) = scounts(pPtr%id) + 1
            pPtr => pPtr%next
         enddo
      endif       

#ifdef _DEBUG_
      if (sum(scounts).ne.npart) then
         STOP 'Error: sum(scounts) .ne. npart!'
      endif
#endif

      !PRINT*,scounts
      !stop
 
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

      ! start at the head of the linked-list
!      if (associated(head)) pPtr => head

      ! sanity check
      !if (.not.associated(head) .and. my_nsend.ne.0) STOP 'WHY?'

      ! SEND

      if (allocated(sbuf)) deallocate(sbuf)
      allocate(sbuf(sum(scounts)))

      do k=1,my_nsend

         m = my_ps(k)

         itag = 1000 + m

#ifdef _DEBUG_
         PRINT*,'rank',id,'sending',scounts(m),'particles to rank',m
#endif

#ifdef _DEBUG_
         if (.not.associated(head)) STOP 'Why is head not associated here?'
#endif

         i = 0
         pPtr => head
         do while (associated(pPtr))

            ! copy linked-list to send buffer
            if (pPtr%id.eq.m) then
               i = i + 1
               sbuf(sdis(k)+i)%tleft = pPtr%tleft
            endif

            pPtr => pPtr%next
         enddo

#ifdef _DEBUG_
         if (i.ne.scounts(m)) STOP 'WHY?'
#endif

         call MPI_SEND(sbuf(sdis(k)+1),scounts(m),MPI_PART_XFT, &
                       m,itag,icomw,ier)

         !do i=1,scounts(m)
         !   PRINT*,sbuf(i)%pid
         !enddo 

      enddo

      ! CLEAR LINKED LIST
      do while (associated(head))
         !!!call remove_link(head) does not work!
         !!!call put_in_hat(head)  does not work!
         pPtr => head
         call remove_link(pPtr) 
         call put_in_hat(pPtr)
      enddo

      ! RECEIVE

      if (allocated(rbuf)) deallocate(rbuf)
      allocate(rbuf(sum(rcounts)))

      do k=1,my_nrecv

         m = my_pr(k)

         itag = 1000 + id

#ifdef _DEBUG_
         PRINT*,'rank',id,'receiving',rcounts(m),'particles from rank',m 
#endif

         call MPI_RECV(rbuf(rdis(k)+1),rcounts(m),MPI_PART_XFT, &
                       m,itag,icomw,mpistat,ier)

      enddo

      ! note: rdis is not always zero
      ! shouldn't need rdis after this point...

      ! here I know I am only sending from 0 to 1
      ! but eventually I will need to break up the
      ! linked-list and copy to the send buffer
      ! in packets, for different ranks 

      !PRINT*,id,'|',size(rbuf)


      ! finally, copy receive buffer to linked-lists
      do i=1,size(rbuf)

         pid = pull_from_hat()
         pPtr => particles(pid)%p
         call insert_link(pPtr)

         pPtr%id = id
         pPtr%tleft = rbuf(i)%tleft      

         !PRINT*,id,'|',rbuf(i)%tleft
      enddo

      ! print in order
      do m=0,nproc-1

         if (id.eq.m) then
            pPtr => head
            do while (associated(pPtr))
               PRINT*,id,'|',pPtr%tleft
               pPtr => pPtr%next
            enddo
         endif
       
         call MPI_BARRIER(icomw,ier)
 
      enddo

      call MPI_Finalize(ier)
    
      end program

