! This program will calculate pref binding coff using equation 56 from DOI: 10.1007/s12013-007-9005-0
program pref_binding_gamma

! intilization
integer :: num, num1 , tot_cro , tot_wat
real, allocatable :: r(:) , cn_cro(:), r1(:) , cn_wat(:)
real , allocatable :: tao(:)


num = 0
num1 = 0

! input of num wat and num cro
tot_cro = 0
tot_wat = 0

write(*,*) "Enter the number of crowders in the system"
read(*,*) tot_cro

write(*,*) "Enter the number of water in the system"
read(*,*) tot_wat


! opening of file
open(unit = 10, file = './cn_crowder.xvg' , status="old", action="read") ! input file crowder
open(unit = 20, file = './cn_water.xvg' , status="old", action="read") ! input file water
open( unit = 50, file = 'output_pref_binding.txt' , action= 'write', status = 'replace') ! outpuut file


! counting number of lines/frames
call count_number_of_lines(10,num)
print*, "Number of lines in your crowder cn file is:", num

call count_number_of_lines(20,num1)
print*, "Number of lines in your water cn file is:", num1

if(num .ne. num1) print*, " Check the input files"
!print*, "You are good to go"


! rewinding to take pointer to initial because we already read the file
rewind(10)
rewind(20)


! allocatable
allocate(r(num) , cn_cro(num), r1(num) , cn_wat(num))
allocate(tao(num))

! reading of cn file
! skipping starting 25 lines
r(:) = 0.0
r1(:) = 0.0
cn_cro(:) = 0.0
cn_wat(:) = 0.0


do i = 1,25
       read(10,*)
       read(20,*) 
end do

do i = 1, (num-25)
        read(10,*) r(i) , cn_cro(i)
        read(20,*) r1(i) , cn_wat(i)
end do        

print*, "Input file is read, now calculating pref_binding"

! calculation of pref binding
tao(:) = 0.0

do j = 1, (num-25)
        tao(j) = cn_cro(j) - ( ((tot_cro - cn_cro(j))/(tot_wat - cn_wat(j))) * cn_wat(j))
        write(50,*) r(j) , tao(j)
end do

!print*, "Everything looks fine"
print*, "Output file is output_pref_binding.txt"



end program pref_binding_gamma

subroutine count_number_of_lines(unit,num)
integer :: unit, num, a
num = -1 ! if 0 it will assign +1 value to EOF also
a = 0
rewind(unit)
do while(a==0)
read(unit,*,iostat=a)
num = num + 1
end do
end subroutine count_number_of_lines

