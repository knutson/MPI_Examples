      program main

      use mpi

      implicit none
      integer, parameter :: NUM_PACKS = 20
      integer :: i,icomw,id,ier,nproc,m,itag,mpistat(MPI_STATUS_SIZE)
      integer, allocatable, dimension(:) :: ir,q
      logical :: done

      call MPI_Init(ier)

      icomw = MPI_COMM_WORLD

      call MPI_Comm_rank(icomw,id,ier)
      call MPI_Comm_size(icomw,nproc,ier)

      if (nproc<2) STOP

      allocate(ir(0:2*nproc-1))
      ir = MPI_REQUEST_NULL

      allocate(q(1))
      q = 0

      if (id.eq.0) then ! send 10 packets to rank 1 

      done = .true.   ! enter first time
      m = 1           ! only send to rank 1
      itag = 1000 + m ! constant + id of receiver
      i = 0

      do while (i<NUM_PACKS)
          
          if (done) then
             i = i + 1
             call MPI_ISEND(i,1,MPI_INTEGER, &
                            m,itag,icomw,ir(m),ier)
             done = .false.
          endif

          call MPI_TEST(ir(m),done,mpistat,ier)
          if (done) then
             PRINT*,'sent packet #',i
             ! free request?
          else
             PRINT*,'communication of packet-to-peer not yet completed'
          endif

      enddo
 
      endif
 
      !call MPI_BARRIER(icomw,ier) ! testing 

      if (id.eq.1) then ! receive 10 packets from rank 0

      done = .true.    ! enter first time 
      m = 0            ! only receive from 0
      itag = 1000 + id ! constant + id of receiver
      i = 0

      do while (i<NUM_PACKS)
          
         if (done) then
            i = i + 1
            call MPI_IRECV(q,1,MPI_INTEGER, &
                           m,itag,icomw,ir(nproc+m),ier)
            done = .false.
         endif

         call MPI_TEST(ir(nproc+m),done,mpistat,ier)

         if (done) then
            ! free request?
            PRINT*,'received packet #',i
            if (i.ne.q(1)) STOP ! sanity check
         else
            PRINT*,'communication of packet-from-peer not yet completed'
         endif

      enddo
     
      endif

      call MPI_WAITALL(2*nproc,ir,MPI_STATUSES_IGNORE,ier)

      deallocate(ir,q)
 
      call MPI_Finalize(ier)
    
      end program
