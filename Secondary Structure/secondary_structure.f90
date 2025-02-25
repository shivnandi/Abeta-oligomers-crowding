! This program will evaluate the data obtained after calculating secondary structure from VMD
program calcalate_ss
character(len=50) :: res_num1, x1, char11
character(len=50) :: res_num2, x2, char22
integer :: n_frames, n_res, frame_num1,frame_num2
character(len=50) , allocatable :: struct1(:),struct2(:,:)
real, allocatable ::  r_coil1(:,:),r_extended1(:,:),r_rurn1(:,:),r_alpha_helix1(:,:)
character(len=2) :: typeA, typeB, typeC , typeD , typeE , typeF , typeG
real ::r1,r4
character(len=50):: r2,r3,a(100)


! input values
write(*,*) " Number of frames are: "
read(*,*) n_frames

write(*,*) " Number of residues are: "
read(*,*) n_res

!allocation
allocate(struct1(n_res),struct2(n_frames,n_res))
allocate(r_coil1(n_frames,n_frames),r_extended1(n_frames,n_frames),r_rurn1(n_frames,n_frames))
allocate(r_alpha_helix1(n_frames,n_frames))



!opening of file

open(unit=10, file= "ss.txt" , status="old", action="read")

!!! reading of file
do i = 1,187
  read(10,*)
end do


do j = 1, n_frames
   do i = 1,4
!! protein chain 1
   do k = 1,n_res
   read(10,*) res_num1, x1, char11, frame_num1 , struct1(k)
  !  call calculate_ss_ratio(struct1, n_res, outextended, outcoil, outturn)
   end do
!print*, n_res   
   call calculate_ss_ratio(struct1, n_res, outextended, outcoil, outturn, outalpha_helix)
 
!    write(*,*) outalpha_helix

    r_coil1(i,j) = outcoil
    r_extended1(i,j) = outextended
    r_rurn1(i,j) = outturn
    r_alpha_helix1(i,j) = outalpha_helix

   write(10+i,*) j, r_coil1(i,j)
   write(20+i,*) j, r_rurn1(i,j)
   write(30+i,*) j, r_extended1(i,j)
   write(40+i,*) j, r_alpha_helix1(i,j)
!   write(*,*) "pritein ends"
   end do
!   write(*,*) "frame ends"
   do m = 1 , 150
   read(10,*)
   end do
!   do n = 1,5
!   read(10,*) r1,r2,r3,r4, a(n)
!   write(*,*) r1,r2,r3,r4,a(n)
!   end do
!   write(*,*)"frame ends"

end do
write(*,*) "Program executed successfully"
end program calcalate_ss


subroutine calculate_ss_ratio(struct1, n_res, outextended, outcoil, outturn, outalpha_helix )
  character(len=50) :: struct1(n_res)
  integer :: n_res
  real :: outextended, outcoil, outturn, outalpha_helix
  integer :: j, l
  real :: turn1, extended1, coil1
  real:: alpha_helix1, bridge1,helix3101,pihelix1
  character(len=2) :: typeA, typeB, typeC , typeD , typeE , typeF , typeG
  
! defining the variables of secondary structure used by STRIDE software
typeA = "T" ! turn
typeB = "H" !alpha helix
typeC = "E" ! extended conformation
typeD = "B" ! isolated bridge
typeE = "G" !3-10 Helix
typeF = "I" !pi-helix
typeG = "C" !coil

  turn1 = 0.0
  extended1 = 0.0
  coil1 = 0.0
  alpha_helix1 = 0.0

  do j = 1, n_res
      if (struct1(j) .eq. typeA) then
        turn1 = turn1 + 1
      elseif (struct1(j) .eq. typeB) then
        alpha_helix1 = alpha_helix1 + 1
      elseif (struct1(j) .eq. typeC) then
        extended1 = extended1 + 1
      elseif (struct1(j) .eq. typeD) then
        bridge1 = bridge1 + 1
      elseif (struct1(j) .eq. typeE) then
        helix3101 = helix3101 + 1
      elseif (struct1(j) .eq. typeF) then
        pihelix1 = pihelix1 + 1
      elseif (struct1(j) .eq. typeG) then
        coil1 = coil1 + 1
      endif
  end do
!  print*, alpha_helix1 , nres
  outalpha_helix = alpha_helix1 / n_res
  outturn = turn1 / n_res
  outextended = extended1 / n_res
  outcoil = coil1 / n_res
!  print*, outalpha_helix
end subroutine   calculate_ss_ratio
