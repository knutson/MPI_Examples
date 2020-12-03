
! TODO: add a routine to assign new particle GIDs
!       add particle diameter, temperature, density?
!       make doubly linked list
!       make particles an array (sized to number of cells)

      module llist

      implicit none
      type, public :: part_t
         integer :: pid = 0       ! (global) particle index
         integer :: i = 0         ! local cell index that contains the particle
         real*8 :: xyz(3) = 0.0d0 ! particle position
         real*8 :: uvw(3) = 0.0d0 ! particle velocity
         real*8 :: tleft  = 0.0d0 ! remaining time until particle has moved dt 
         type(part_t), pointer :: next=>NULL()
         !type(part_t), pointer :: prev=>NULL()
      end type

      type, public :: partlist_t
         type(part_t), pointer :: p=>null()
         integer :: count = 0
      end type

      !global
      type(part_t), pointer :: head=>null()
      !type(part_t), pointer :: tail=>null()

      type(partlist_t) :: particles
      !type(partlist_t), allocatable :: part_lists(:)

      contains

      subroutine insert_item(item)

      implicit none
      type(part_t), intent(IN) :: item
      type(part_t), pointer :: pPtr

      if (.not.associated(head)) then ! list is empty, make item the head
         allocate(particles%p)
         head => particles%p
         head%pid = item%pid
      else

         pPtr => head

         ! check head
         if (pPtr%pid.eq.item%pid) then
            STOP 'Error: attempted to insert item that is already in the list!'
         endif

         do while (associated(pPtr%next))

            pPtr => pPtr%next

            if (pPtr%pid.eq.item%pid) then
               STOP 'Error: attempted to insert item that is already in the list!'
            endif

         enddo

         allocate(pPtr%next)
         pPtr => pPtr%next
         pPtr%pid = item%pid

      endif

      particles%count = particles%count + 1

      end subroutine

      function get_wtime() result(wtime)
      implicit none
      real*8 :: wtime
      integer(kind=4) :: h,m,s,ms,values(8)
      call date_and_time( values = values )
      h = values(5)
      m = values(6)
      s = values(7)
      ms = values(8)
      wtime = dble( 3600*h + 60*m + s ) + dble(ms)/1000.0d0
      end function 

      subroutine print_list
      implicit none
      type(part_t), pointer :: pPtr 

      if (.not.associated(head)) then
         PRINT*,'List is empty!'
      else
         pPtr => head
         do while (associated(pPtr))
            write(6,*) pPtr%pid
            pPtr => pPtr%next
         enddo
         !note: pPtr is not associated at this point
      endif
      !PRINT*

      end subroutine

      subroutine delete_item(pid)

      implicit none
      integer, intent(in) :: pid
      type(part_t), pointer :: pPtr,dPtr
      logical :: found

      if (.not.associated(head)) then
         PRINT*,'WARNING: cannot delete an item from an empty list!'
         !STOP
      else

         if (head%pid.eq.pid) then ! we are deleting the first item

            dPtr => head
            head => head%next ! null if only 1 item in the list
            deallocate(dPtr) ! CRUCIAL
            particles%count = particles%count - 1

         else

            found = .false.
            pPtr => head
            do while (associated(pPtr%next))         
               if (pPtr%next%pid.eq.pid) then
                  found = .true.
                  exit
               else
                  pPtr => pPtr%next
               endif
            enddo
   
            if (found) then
               ! note: the item to be deleted is pPtr%next
               dPtr => pPtr%next
               if (.not.associated(pPtr%next%next)) then ! we are deleting the last item
                  pPtr%next => null()
                  !nullify( pPtr%next )
               else
                  pPtr%next => pPtr%next%next
               endif
               deallocate(dPtr) ! CRUCIAL
               particles%count = particles%count - 1

            else
               PRINT*,'WARNING: Could not find item to delete!'
            endif         

         endif

      endif

      end subroutine

      end module

