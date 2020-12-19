
      program main

      use llist

      implicit none

      integer, parameter :: id = 0

      integer :: i,pid
      type(part_t), pointer :: pPtr
      type(part_t), pointer :: pPtr_next

      ! initialize the hat (this is where PIDs are set)
      call init_hat()

      do i=1,11
         pid = pull_from_hat()
         pPtr => particles(pid)%p
         call insert_link(pPtr)

         if (mod(i,2).ne.0) then ! i is odd
            pPtr%id = 1
         else                    ! i is even
            pPtr%id = 0
         endif
      enddo
 
      pPtr => head
      do while(associated(pPtr))

         if (pPtr%id.ne.id) then
            pPtr_next => pPtr%next
            call remove_link(pPtr)
            call put_in_hat(pPtr)
            pPtr => pPtr_next
         else
            pPtr => pPtr%next 
         endif

      enddo

      ! print
      pPtr => head
      do while(associated(pPtr))
         PRINT*,pPtr%pid
         pPtr => pPtr%next 
      enddo
    
      end program

