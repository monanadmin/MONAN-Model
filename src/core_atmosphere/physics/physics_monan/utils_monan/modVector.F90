module ModVector

   implicit none

   private

   character(len=*), parameter   :: p_source_name = 'modVector.F90'
   
   
   ! Vector inside data type
   type :: data_t
      integer :: value
   end type data_t
   ! Vector data type
   
   type :: vector_t
      private
      type(data_t), pointer, dimension(:) :: vector
      integer :: num_elements
      
   end type vector_t

   public :: vector_t
   public :: init
   public :: free_memory
   public :: get_size
   public :: get_num_elements
   public :: insert
   public :: vector_put
   public :: get_data
   public :: get_data_value
   public :: remove
   public :: print_all
   public :: insert_range
   public :: insert_unique


contains

  ! Initialize vector
   subroutine init(instance, vec_max_size)
      implicit none 
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: vec_max_size

      instance%num_elements = 0
      allocate(instance%vector(vec_max_size))
   end subroutine init


  ! Insert a range of values in vector
   subroutine insert_range(instance, start_val, end_val)
      implicit none
      type(vector_t), intent(inout) :: instance
      
      integer, intent(in)           :: start_val
      integer, intent(in)           :: end_val
      integer                       :: i, vector_ind
      character(len=*), parameter   :: p_procedure_name = 'insert_range' 
      
      if (end_val < start_val) then
         print*, p_procedure_name, "Erro: end_val", end_val, " menor que start_val", start_val 
         stop
       endif
            
      vector_ind=0      
      do i = start_val, end_val
         vector_ind=vector_ind+1
         instance%vector(vector_ind)%value=i
      enddo
      instance%num_elements = vector_ind
  
   end subroutine insert_range


  ! Free the entire list and all data, beginning at SELF
   subroutine free_memory(instance)
      implicit none
      type(vector_t), intent(inout) :: instance
      if(associated(instance%vector)) deallocate(instance%vector)
      instance%num_elements = 0
   end subroutine free_memory


   ! get vector max size (by init)
   function get_size(instance) result(size_vec)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer :: size_vec
      if(associated(instance%vector)) then
         size_vec = size(instance%vector)
      else
         size_vec = 0
      endif
   end function get_size


   ! get vector actual num of elements
   function get_num_elements(instance) result(num_elements)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer :: num_elements
      num_elements = instance%num_elements
   end function get_num_elements


   ! Insert a value at end of the vector
   subroutine insert(instance, value_param)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: value_param

      instance%num_elements = instance%num_elements + 1
      instance%vector(instance%num_elements)%value = value_param
   end subroutine insert


   ! Insert a value at end of the vector, if not exists
   function insert_unique(instance, value_param) result(is_inserted)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: value_param
      logical :: is_inserted
      integer :: index

      is_inserted = .false.

      index = get_index_of_value(instance, value_param)
      if(index == -1) then 
         instance%num_elements = instance%num_elements + 1
         instance%vector(instance%num_elements)%value = value_param
         is_inserted = .true.
      endif
      
   end function insert_unique


   ! Store the encoded data the index_data position
   subroutine vector_put(instance, value_param, data_index)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: value_param
      integer, intent(in) :: data_index

      instance%vector(data_index)%value = value_param
   end subroutine vector_put


   ! Return the DATA stored in data_index
   function get_data(instance, data_index) result(data)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: data_index
      type(data_t) :: data
      data = instance%vector(data_index)
   end function get_data


   function get_data_value(instance, data_index) result(data_value)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: data_index
      integer :: data_value
      data_value = instance%vector(data_index)%value
   end function get_data_value


   function get_index_of_value(instance, data_value) result(index)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: data_value
      integer :: index, i

      index = -1
      do i = 1, instance%num_elements
         if (instance%vector(i)%value == data_value) then
            index = i
         endif
      enddo 

   end function
   
   ! print all elements of vector
   subroutine print_all(instance, views) 
      implicit none
      type(vector_t), intent(inout) :: instance
      integer :: data_index
      integer, intent(in) :: views
      
      write(*, '(A8)', advance='NO') 'vector = ('
      do data_index = 1, min(instance%num_elements, views)
         write(*,'(i8, "," )',advance='NO')  instance%vector(data_index)%value
      end do
      if (instance%num_elements > views) then
         ! print last elements 
         write(*,'(A5)',advance='NO') ' ... '
         do data_index = max(instance%num_elements - views, views +1 ), instance%num_elements
         write(*,'(i8, "," )',advance='NO')  instance%vector(data_index)%value
         end do
      endif
      write(*, '(A2)', advance='YES') ' )'

   end subroutine print_all


   ! remove element containing data from parameter
   function remove(instance, value_param) result(is_removed)
      implicit none
      type(vector_t), intent(inout) :: instance
      integer, intent(in) :: value_param
      integer :: index
      logical :: is_removed

      is_removed = .false.
      if (instance%num_elements == 0) then
         return
      endif

      index = get_index_of_value(instance, value_param)
      if(index == -1) return

      instance%vector(index:instance%num_elements-1) = instance%vector(index+1:instance%num_elements)
      instance%num_elements = instance%num_elements -1
      is_removed = .true.

   end function remove

end module ModVector