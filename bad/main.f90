      program main

      use mpi

      implicit none

      integer :: icomw,id,ier,nproc,m,itag
      integer, allocatable, dimension(:) :: ir,q

      ! all processes send their id
      ! each process receives ordered id list in buffer q
      ! for example, if nproc = 4, q = 0 1 2 3 on every process

      logical :: flag,done

      integer :: mpistat(MPI_STATUS_SIZE)

      call MPI_Init(ier)

      icomw = MPI_COMM_WORLD

      call MPI_Comm_rank(icomw,id,ier)
      call MPI_Comm_size(icomw,nproc,ier)

      allocate( ir(0:2*nproc-1), q(nproc) )

      ir = MPI_REQUEST_NULL

      q = 0

      do m = 0,nproc-1
         itag = 1000 + m ! constant + id of receiver
         call MPI_ISEND(id,1,MPI_INTEGER, &
              m,itag,icomw,ir(m),ier)
      enddo

      do m = 0,nproc-1
         itag = 1000 + id ! constant + id of receiver
         call MPI_IRECV(q(m+1),1,MPI_INTEGER, &
              m,itag,icomw,ir(nproc+m),ier)
      enddo

      done = .false.
      do while (.not.done)

         done = .true. ! unless any messages have not been received 

         do m=0,nproc-1

            ! check receives 
            !call MPI_TEST(ir(nproc+m),flag,mpistat,ier) ! deallocated request 
            call MPI_Request_get_status(ir(nproc+m),flag,mpistat,ier) ! non-destructive test
      
            if (flag) then
               ! mpistat is only significant if message has been received (flag=.true.)?
               !PRINT*,id,m,mpistat(MPI_SOURCE)
               if (m.ne.mpistat(MPI_SOURCE)) then
                  PRINT*,'why did this happen?',m,mpistat(MPI_SOURCE)
                  !STOP 
               endif
            else
               done = .false.
               !PRINT*,'not done yet'
            endif

         enddo

      enddo
      PRINT*,'rank',id,'has received all messages'

      call MPI_WAITALL(2*nproc,ir,MPI_STATUSES_IGNORE,ier)

      deallocate(ir,q)
 
      call MPI_Finalize(ier)
    
      end program
