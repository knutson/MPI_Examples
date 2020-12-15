
      ! we are assuming calls to pull_from_hat is always followed by insert_link
      ! and remove_link is followed by put_in_hat
      ! so why not combine? 

      module llist

      implicit none
      type, public :: part_t
         integer :: pid     = -1      
         integer :: i       = -1          ! local index of cell that contains the particle
         integer :: id      = -1          ! rank which owns cell that contains the particle
         real*8 :: xyz(3)   = huge(1.0d0) ! particle position
         real*8 :: uvw(3)   = huge(1.0d0) ! current particle velocity
         real*8 :: uvw_n(3) = huge(1.0d0) ! previous particle velocity
         real*8 :: tleft    = huge(1.0d0) ! remaining time until particle has moved pdt (particle dt) 
         type(part_t), pointer :: next=>null()
         type(part_t), pointer :: prev=>null()
      end type

      type, public :: partlist_t
         type(part_t), pointer :: p=>null()
      end type

      ! global

      type(part_t), pointer :: head=>null()
      type(part_t), pointer :: tail=>null()

      type(part_t), pointer :: phat=>null()

      type(partlist_t), allocatable :: particles(:)

      integer :: npart ! local particle count

      contains

      ! create (expand) hat
      subroutine init_hat()
 
      implicit none
      integer :: i

      npart = 0

      allocate(particles(10000))
      do i=10000,1,-1 
         allocate(particles(i)%p)
#ifdef _DEBUG_
         call clear_part(particles(i)%p) ! Is this necessary? If so, when? During expansion?
#endif
         particles(i)%p%pid = i
         if (associated(phat)) particles(i)%p%next => phat ! what does this do?
         phat => particles(i)%p
      enddo

      end subroutine

      function pull_from_hat() result(pid)
      implicit none
      integer :: pid
      if (.not.associated(phat)) then ! need to expand hat, give error for now
         STOP 'Error: hat is not big enough'
      endif         
      pid = phat%pid
      phat => phat%next
      npart = npart + 1
      nullify(particles(pid)%p%next) ! shouldn't this already be null?
      end function

      subroutine put_in_hat(pPtr)
      implicit none
      type(part_t), pointer :: pPtr
      call clear_part(pPtr)
      if (associated(phat)) pPtr%next => phat
      phat => pPtr 
      npart = npart - 1 
      end subroutine

      ! insert link at the end of the linked-list (push back)
      subroutine insert_link(pPtr)

      implicit none
      type(part_t), pointer :: pPtr

      if (.not.associated(head)) then
         head => pPtr
#ifdef _DEBUG_
         if (associated(pPtr%next)) STOP 'How did this happen? (1)'
         if (associated(pPtr%prev)) STOP 'How did this happen? (2)'
#endif
      endif
      if (associated(tail)) then
         ! note: here the old tail is becoming the second to last item
         tail%next => pPtr
         pPtr%prev => tail
      endif
      ! update the tail
      tail => pPtr
     
      end subroutine

      subroutine remove_link(pPtr)
      implicit none
      type(part_t), pointer, intent(inout) :: pPtr
      logical :: lprev,lnext

      lprev = associated(pPtr%prev)
      lnext = associated(pPtr%next)

      if (lprev.and.lnext) then ! neither head nor tail
         pPtr%prev%next => pPtr%next
         pPtr%next%prev => pPtr%prev
      else if (lnext) then ! remove the head
         head => pPtr%next
         nullify(head%prev)
      else if (lprev) then ! remove the tail
         tail => pPtr%prev
         nullify(tail%next)
      else ! remove the only link
         head => null()
         tail => null()
      endif

      end subroutine

      subroutine clear_part(pPtr)
      implicit none
      type(part_t), pointer :: pPtr
      
      !NOTE: DO NOT CLEAR PID!
      pPtr%i        = -1          ! local index of cell that contains the particle
      pPtr%id       = -1          ! rank which owns cell that contains the particle
      pPtr%xyz(3)   = huge(1.0d0) ! particle position
      pPtr%uvw(3)   = huge(1.0d0) ! current particle velocity
      pPtr%uvw_n(3) = huge(1.0d0) ! previous particle velocity
      pPtr%tleft    = huge(1.0d0) ! remaining time until particle has moved pdt (particle dt) 
    
      nullify(pPtr%prev,pPtr%next)
      end subroutine 

      end module

